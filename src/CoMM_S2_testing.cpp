//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "math.h"
#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include <stdio.h>
#include <bitset>
#include <iostream>
#include <fstream>
#include <sstream>
#include "SCLMM.hpp"
#include <boost/algorithm/string.hpp>
#include <R.h>
#include "readPlink2.hpp"
#include "StandardizeData.hpp"
#include "mt_paral_testing_job2.hpp"
#include "readExprFile.hpp"


#define MAX_LEN 20

using namespace std;
using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH, bigmemory)]]
// [[Rcpp::plugins(cpp11)]]


List read_GWAS(string filename, int P){
  std::ifstream ifs(filename.c_str());
  std::string line;
  vector<string> fields;
  CharacterVector snp(P - 1);
  IntegerVector A1(P - 1);
  IntegerVector A2(P - 1);
  uvec chr = zeros<uvec>(P - 1);
  uvec BP = zeros<uvec>(P - 1);
  mat GWAS_sum = zeros(P - 1, 2); // transpcriptome summary data
  int i = 0;
  string str;
  while (std::getline(ifs, line)) // read one line from ifs
  {
    std::istringstream iss(line); // access line as a stream
    boost::split(fields, line, boost::is_any_of(" \t *"));
    
    if (i > 0){
      snp[i-1] = fields[0];
      chr[i-1] = atof(fields[1].c_str());
      BP[i-1] = atof(fields[2].c_str());
      str = fields[3];
      char* chr1 = strdup(str.c_str());
      A1[i-1] = (int)(*chr1);
      str = fields[4];
      char* chr2 = strdup(str.c_str());
      A2[i-1] = (int)(*chr2);
      GWAS_sum(i-1, 0) = atof(fields[5].c_str());
      GWAS_sum(i-1, 1) = atof(fields[6].c_str());
    }
    i++;
  }
  ifs.close();
  
  List output = List::create(
    Rcpp::Named("GWAS_snps") = snp,
    Rcpp::Named("GWAS_chr") = chr,
    Rcpp::Named("BP") = BP,
    Rcpp::Named("A1") = A1,
    Rcpp::Named("A2") = A2,
    Rcpp::Named("GWAS_sum") = GWAS_sum);
  return output;
}


Rcpp::List ReadSNPinfo(std::string stringname, IntegerVector A1, IntegerVector A2, CharacterVector rsname,
                       IntegerVector chr, IntegerVector bp, NumericVector morgan, int N)
{
  FILE *stream;
  int ch, b;
  double mor;
  char s[MAX_LEN + 1], efa, nefa;
  
  stream = fopen(stringname.c_str(), "r");
  clock_t t1 = clock();
  //int i = 0;
  /* Put in various data. */
  for (int i = 0; i < N; i++){
    if (i % 100000 == 0 && i != 0){
      cout << i << "-th SNP" << ",";
      cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
    }
    
    fscanf(stream, "%i %s %lf %i %c %c", &ch, &s[0], &mor, &b, &efa, &nefa);
    
    chr(i) = ch;
    rsname(i) = s;
    morgan(i) = mor;
    bp(i) = b;
    A1(i) = (int)efa;
    A2(i) = (int)nefa;
    
  }
  List output = List::create(Rcpp::Named("chr") = chr,
                             Rcpp::Named("rsname") = rsname,
                             Rcpp::Named("morgan") = morgan,
                             Rcpp::Named("bp") = bp,
                             Rcpp::Named("A1") = A1,
                             Rcpp::Named("A2") = A2);
  return output;
}

void T2AG2C(uvec& A31_){
  
  uvec idx;
  idx = find(A31_ == 84);
  uvec idxrepl1(idx.n_elem);
  idxrepl1.fill(65);
  A31_.elem(idx) = idxrepl1;
  
  idx = find(A31_ == 71);
  uvec idxrepl2(idx.n_elem);
  idxrepl2.fill(67);
  A31_.elem(idx) = idxrepl2;
  
}


Rcpp::List getColNum_Header2(std::string filename, char delimiter){
	//cout << "Count columns in the file: ";

	std::ifstream myfile(filename.c_str());
	std::string line, temp;

	if (myfile.is_open()){
		getline(myfile, line);
	}

	stringstream ss(line); // Turn the string into a stream.
	string tok;
	CharacterVector header;

	while (getline(ss, tok, delimiter)) {
		header.push_back(tok);
	}

	//cout << "line 1" << line;
	std::istringstream iss(line);
	int columns = 0;
	do{
		std::string sub;
		iss >> sub;
		if (sub.length())
			++columns;
	} while (iss);

	//cout << columns << endl;

	List output = List::create(//Rcpp::Named("NumGene") = NumGene,
		Rcpp::Named("columns") = columns,
		Rcpp::Named("header") = header);

	return output;
}


CharacterVector charv_subset2_(CharacterVector x, uvec idx){
	CharacterVector v(idx.n_elem);
	for (unsigned int i = 0; i < idx.n_elem; i++){
		v(i) = x(idx(i));
	}
	return v;
}


void ReadPlinkFamFile2(std::string stringname, CharacterVector FID, CharacterVector IID, IntegerVector sex,
	NumericVector pheno, int N)
{
	FILE *stream;
	int gender;
	double phn;
	char fid[MAX_LEN + 1], iid[MAX_LEN + 1], tmp1[MAX_LEN + 1], tmp2[MAX_LEN + 1];

	stream = fopen(stringname.c_str(), "r");

	for (int i = 0; i < N; i++){
		fscanf(stream, "%s %s %s %s %i %lf", &fid, &iid, &tmp1, &tmp2, &gender, &phn);

		FID(i) = fid;
		IID(i) = iid;
		sex(i) = gender;
		pheno(i) = phn;

	}
}


List cv_glmnet(mat x, vec y, double alpha){

	// calling cv.glmnet()
	Function f("cv.glmnet");
	return f(Named("x") = x, Named("y") = y, Named("alpha") = alpha);

}

List glmnet(mat x, vec y, double lambda, double alpha){

	// calling glmnet()
	Function f("glmnet");
	return f(Named("x") = x, Named("y") = y, Named("lambda") = lambda, Named("alpha") = alpha);

}



// [[Rcpp::export]]
List CoMM_S2_paral_testing(std::string stringname1, std::string stringname2, std::string stringname3, std::string stringname4, std::string stringname5,
	int bw = 500000, double lam = 0.95, const int coreNum = 1){
	std::string method = "pdsce";
	// Genotype for eQTL data: stringname1; 
	// Summary data for GWAS data; stringname2;
	// Individual data for reference panel; stringname3;
	// expression file for eQTL data; stringname4;
	// covariate file stringname5;

	//-------------------------------------------------------------------------------------------//
	cout << "testing in parallel " << endl;
	cout << "begin to read eQTL data" << endl;

	string famfile1 = stringname1;
	famfile1 += ".fam";
	int N1 = getLineNum(famfile1);
	string bimfile1 = stringname1;
	bimfile1 += ".bim";
	long int P1 = getLineNum(bimfile1);

	IntegerVector A11(P1), A12(P1);
	CharacterVector rsname1(P1);
	IntegerVector chr1(P1), bp1(P1);
	NumericVector morgan1(P1);

	//read from SNPinfo file (pass as pointer)
	cout << "## Start loading SNP info:" << endl;
	clock_t t1 = clock();
	ReadSNPinfo(bimfile1, A11, A12, rsname1, chr1, bp1, morgan1, P1);
	cout << "Finish loading SNP info in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;

	// load IID in file 1 (fam file 1)
	cout << "## Start loading fam files 1 " << endl;
	IntegerVector tmp1(N1);
	CharacterVector FID_1(N1), IID_1(N1);
	NumericVector tmp2(N1);
	ReadPlinkFamFile(famfile1, FID_1, IID_1, tmp1, tmp2, N1);
	cout << "Finish loading fam files 1 in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;

	ObjXY obj1_XY = ReadDataFromFile(stringname1);
	Mat<unsigned>* X1_point = &obj1_XY.X;


	//-------------------------------------------------------------------------------------------//
	cout << "begin to read GWAS data" << endl;

	int P2 = getLineNum(stringname2);
	List GWAS = read_GWAS(stringname2, P2);

	IntegerVector A21(P2), A22(P2);
	CharacterVector rsname2(P2);
	IntegerVector chr2(P2), bp2(P2);

	A21 = GWAS["A1"];
	A22 = GWAS["A2"];
	rsname2 = GWAS["GWAS_snps"];
	chr2 = GWAS["GWAS_chr"];
	bp2 = GWAS["BP"];
	mat GWAS_sum = GWAS["GWAS_sum"];


	//-------------------------------------------------------------------------------------------//
	cout << "begin to read reference panel" << endl;

	string famfile3 = stringname3;
	famfile3 += ".fam";
	int N3 = getLineNum(famfile3);
	string bimfile3 = stringname3;
	bimfile3 += ".bim";
	long int P3 = getLineNum(bimfile3);

	IntegerVector A31(P3), A32(P3);
	CharacterVector rsname3(P3);
	IntegerVector chr3(P3), bp3(P3);
	NumericVector morgan(P3);

	//read from SNPinfo file (pass as pointer)
	cout << endl;
	cout << "Start loading SNP info:" << endl;
	clock_t t3 = clock();
	ReadSNPinfo(bimfile3, A31, A32, rsname3, chr3, bp3, morgan, P3);
	cout << "Finish loading SNP info in " << (clock() - t3)*1.0 / CLOCKS_PER_SEC << " sec." << endl;

	ObjXY obj3_XY = ReadDataFromFile(stringname3);
	Mat<unsigned>* X3_point = &obj3_XY.X;


	//-------------------------------------------------------------------------------------------//
	cout << "## Start loading expression file. " << endl;
	int Ngene = getLineNum(stringname4) - 1; // header is first line;
	char delimiter = '\t';
	List tmp = getColNum_Header2(stringname4, delimiter);

	int Nindiv = tmp["columns"];
	CharacterVector header = tmp["header"];

	// IID from header of expr file
	uvec idxtmp = linspace<uvec>(6, Nindiv - 1, Nindiv - 6);
	CharacterVector IID_expr = charv_subset2_(header, idxtmp);

	// match file 1 with expression file
	CharacterVector indiv_inter = intersect(IID_1, IID_expr);
	IntegerVector idxint = match(IID_1, indiv_inter); // put the IID in the order of file 1
	CharacterVector indiv_4use = IID_1[Rcpp::is_na(idxint) == false];
	IntegerVector idxin1 = match(indiv_4use, IID_1) - 1;  //index for individual in file 1
	IntegerVector idxin2 = match(indiv_4use, IID_expr) - 1;  //index for individual in expr file

	uvec idxin1_ = as<uvec>(idxin1);
	uvec idxin2_ = as<uvec>(idxin2);

	//read expr file
	List expr_tmp = getExprFile2(stringname4, delimiter, Nindiv, Ngene + 1, 6);
	CharacterVector genetype1 = expr_tmp["genetype1"];
	CharacterVector genetype2 = expr_tmp["genetype2"];
	CharacterVector targetID = expr_tmp["targetID"];
	mat expr = expr_tmp["expr"];
	vec lower = expr_tmp["lower"];
	vec upper = expr_tmp["upper"];
	vec chr_expr = expr_tmp["chr_expr"];

	// subset overlapped individuals
	mat* exprTmp;
	mat expr_used;
	exprTmp = &expr;
	expr_used = exprTmp->cols(as<uvec>(idxin2));

	//-------------------------------------------------------------------------------------------//
	cout << "matching SNPs and correcting the direction of minor allele" << endl;

	CharacterVector rs_12 = intersect(rsname1, rsname2);
	CharacterVector rs_123 = intersect(rs_12, rsname3);

	IntegerVector idxin = match(rsname1, rs_123);
	CharacterVector rsname_4use = rsname1[Rcpp::is_na(idxin) == false]; // rsname_4use in the order of eQTL.
	IntegerVector bp1_4use = bp1[Rcpp::is_na(idxin) == false]; // bp1 in the order of eQTL
	IntegerVector chr1_4use = chr1[Rcpp::is_na(idxin) == false]; // chr1 in the order of eQTL

	IntegerVector idx_in1 = match(rsname_4use, rsname1);  //index for SNPs in eQTL
	IntegerVector idx_in2 = match(rsname_4use, rsname2);  //index for SNPs in GWAS
	IntegerVector idx_in3 = match(rsname_4use, rsname3);  //index for SNPs in reference panel
	cout << "Size of matched SNPs: " << rsname_4use.size() << endl;

	//----------------------------------------------------------------------------------------------//
	//compare direction
	uvec temp;

	uvec idx1 = as<uvec>(idx_in1) -1;
	temp = as<uvec>(A11);
	uvec A11_ = temp.elem(idx1);
	temp = as<uvec>(A12);
	uvec A12_ = temp.elem(idx1);

	uvec idx2 = as<uvec>(idx_in2) -1;
	temp = as<uvec>(A21);
	uvec A21_ = temp.elem(idx2);
	temp = as<uvec>(A22);
	uvec A22_ = temp.elem(idx2);

	uvec idx3 = as<uvec>(idx_in3) -1;
	temp = as<uvec>(A31);
	uvec A31_ = temp.elem(idx3);
	temp = as<uvec>(A32);
	uvec A32_ = temp.elem(idx3);

	//----------------------------------------------------------------------------------------------//
	//ascii: (A:65; C:67; G:71; T:84) (a:97;c:99,g:103;t:116)
	//replace lower letter to upper letter in A11_, A12_, A21_ and A22_;

	uvec idx;

	idx = find(A11_ > 85);
	A11_.elem(idx) = A11_.elem(idx) - 32;
	idx = find(A12_ > 85);
	A12_.elem(idx) = A12_.elem(idx) - 32;
	idx = find(A11_ > 85);
	A21_.elem(idx) = A21_.elem(idx) - 32;
	idx = find(A12_ > 85);
	A22_.elem(idx) = A22_.elem(idx) - 32;
	idx = find(A31_ > 85);
	A31_.elem(idx) = A31_.elem(idx) - 32;
	idx = find(A32_ > 85);
	A32_.elem(idx) = A32_.elem(idx) - 32;


	//A11_: replace T with A,replace G with C
	T2AG2C(A11_);

	//A12_: replace T with A,replace G with C
	T2AG2C(A12_);

	//A31_: replace T with A, replace G with C
	T2AG2C(A31_);

	//A32_: replace T with A,replace G with C
	T2AG2C(A32_);

	//A21_: replace T with A,replace G with C
	T2AG2C(A21_);

	//A22_: replace T with A,replace G with C
	T2AG2C(A22_);

	//remove index having the same allele
	uvec idxtmp1, idxtmp2;
	idxtmp1 = find((A31_ + A32_) == (A11_ + A12_));
	idxtmp2 = find((A31_ + A32_) == (A21_ + A22_));

	uvec idx4 = intersect(idxtmp1, idxtmp2);

	uvec A11_r = A11_.elem(idx4), A12_r = A12_.elem(idx4);
	uvec A21_r = A21_.elem(idx4), A22_r = A22_.elem(idx4);
	uvec A31_r = A31_.elem(idx4), A32_r = A32_.elem(idx4);

	IntegerVector idx_tmp = wrap(idx4);
	CharacterVector rsname_4use_r = rsname_4use[idx_tmp];
	IntegerVector bp_4use_r = bp1_4use[idx_tmp];
	IntegerVector chr_4use_r = chr1_4use[idx_tmp];

	cout << "Size of A31_r is " << A31_r.n_elem << endl;
	cout << "Size of A11_r is " << A11_r.n_elem << endl;
	cout << "Size of A21_r is " << A21_r.n_elem << endl;

	IntegerVector idxinFile1_ = match(rsname_4use_r, rsname1) - 1;  //index for SNPs in eQTL
	IntegerVector idxinFile2_ = match(rsname_4use_r, rsname2) - 1;  //index for SNPs in GWAS
	IntegerVector idxinFile3_ = match(rsname_4use_r, rsname3) - 1; //index for SNPs in reference panel

	uvec idxinFile1 = as<uvec>(idxinFile1_);
	uvec idxinFile2 = as<uvec>(idxinFile2_);
	uvec idxinFile3 = as<uvec>(idxinFile3_);

	uvec idx21 = find(A11_r != A21_r);
	cout << "Size of A11_r != A21_r is " << idx21.n_elem << endl;
	uvec temp1 = idx2(idx4);
	GWAS_sum(temp1(idx21), zeros<uvec>(1)) = -GWAS_sum(temp1(idx21), zeros<uvec>(1));

	uvec idx22 = find(A11_r != A31_r);
	cout << "Size of A11_r != A31_r is " << idx22.n_elem << endl;
	uvec temp2 = idx3(idx4);
	X3_point->cols(temp2(idx22)) = 2 - X3_point->cols(temp2(idx22));

	//----------------------------------------------------------------------------------------------//

	mat expr_info = zeros<mat>(lower.n_elem, 3);
	expr_info.col(0) = chr_expr;
	expr_info.col(1) = lower;
	expr_info.col(2) = upper;

	vec bp_4use_r_vec = as<vec>(bp_4use_r);
	vec chr_4use_r_vec = as<vec>(chr_4use_r);

	mat snp_info = zeros<mat>(chr_4use_r.length(), 2);
	snp_info.col(0) = bp_4use_r_vec;
	snp_info.col(1) = chr_4use_r_vec;

	//---------------------------------------------------------------//
	// centering
	int size2 = X1_point->n_cols;
	vec meanX(size2);
	mat X1 = zeros(N1, size2);

	for (int i = 0; i < size2; i++){
		meanX[i] = sum(X1_point->col(i))*1.0 / N1;
		vec v_i = conv_to<vec>::from(X1_point->col(i));
		v_i = v_i - meanX[i];
		X1.col(i) = v_i;
	}

	int size3 = X3_point->n_cols;
	vec meanX3(size3);
	mat X3 = zeros(N3, size3);

	for (int i = 0; i < size3; i++){
		meanX3[i] = sum(X3_point->col(i))*1.0 / N3;
		vec v_i = conv_to<vec>::from(X3_point->col(i));
		v_i = v_i - meanX3[i];
		X3.col(i) = v_i;
	}

	//----------------------------------------------------------------//
	// testing
	t1 = clock();
	cout << "### Precompute the active genes ... ";
	uvec idx_active_gene = zeros<uvec>(0);
	uvec g_tmp(1);
	for (uword g = 0; g < Ngene; g++){

		idx = find(bp_4use_r_vec < upper(g) + bw && bp_4use_r_vec > lower(g) - bw
			&& chr_4use_r_vec == chr_expr(g));
		if (idx.is_empty() == false){
			if (idx.n_elem > 1){
				g_tmp(0) = g;
				idx_active_gene = join_cols(idx_active_gene, g_tmp);
			}
		}
	}
	cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;


	// combine gene info into mat
	uword Ngene_active = idx_active_gene.n_elem;
	mat gene_info = zeros<mat>(Ngene_active, 3);
	gene_info = expr_info.rows(idx_active_gene);
	mat expr_active = expr_used.rows(idx_active_gene);

	CharacterVector gene_type1 = genetype1[as<Rcpp::IntegerVector>(wrap(idx_active_gene))];
	CharacterVector gene_type2 = genetype2[as<Rcpp::IntegerVector>(wrap(idx_active_gene))];

	vec stat = zeros(Ngene_active, 1);
	vec alphag = zeros(Ngene_active, 1);
	uvec NumSNP = zeros<uvec>(Ngene_active, 1);

	cout << "### Start running CoMM_S2 on active genes ... " << endl;;
	t1 = clock();

	mat out_param0 = -99 * ones<mat>(Ngene_active, 7);

	//set parallel structure object
	parGene_SCLMM_IS parObj(gene_info, snp_info, bw, N1, P1, N3, P3,
		idxinFile1, idxinFile2, idxinFile3, idxin1_,
		expr_active, X1, X3, GWAS_sum,
		alphag, stat, NumSNP, Ngene_active, lam, method, out_param0);

	const int n_thread = coreNum;
	std::vector<std::thread> threads(n_thread);
	for (int i_thread = 0; i_thread < n_thread; i_thread++){
		threads[i_thread] = std::thread(&parGene_SCLMM_IS::update_by_thread_SCLMM, &parObj, i_thread);
	}

	for (int i = 0; i < n_thread; i++){
		threads[i].join();
	}

	cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;

	mat out_param = -99 * ones<mat>(Ngene, 7);
	out_param.rows(idx_active_gene) = parObj.out_param;

	List output = List::create(
		Rcpp::Named("expr") = expr,  // expression
		Rcpp::Named("expr_active") = expr_active, //expression for active
		Rcpp::Named("rsname_4use_r") = rsname_4use_r,
		Rcpp::Named("snp_info_4use_r") = snp_info,
		Rcpp::Named("targetID") = targetID,
		Rcpp::Named("genetype1") = genetype1,
		Rcpp::Named("genetype2") = genetype2, // all gene name
		Rcpp::Named("gene_type1") = gene_type1,
		Rcpp::Named("gene_type2") = gene_type2, // active gene name
		Rcpp::Named("expr_info") = expr_info, // gene information for all genes: chr lower upper
		Rcpp::Named("gene_info") = gene_info, // gene information for active genes: chr lower upper
		Rcpp::Named("stat") = parObj.stat,
		Rcpp::Named("alphag") = parObj.alphag,
		Rcpp::Named("NumSNP") = parObj.NumSNP,
		Rcpp::Named("out_param") = out_param,
		Rcpp::Named("out_param0") = parObj.out_param);
	return output;

}






























// [[Rcpp::export]]
List CoMM_S2_testing(std::string stringname1, std::string stringname2, std::string stringname3, std::string stringname4, std::string stringname5,
	int bw = 500000, double lam = 0.95, double alpha = 0){

	// Genotype for eQTL data: stringname1; 
	// Summary data for GWAS data; stringname2;
	// Individual data for reference panel; stringname3;
	// expression file for eQTL data; stringname4;
	// covariate file stringname5;

	//-------------------------------------------------------------------------------------------//
	cout << "testing in sequence " << endl;
	cout << "begin to read eQTL data" << endl;

	string famfile1 = stringname1;
	famfile1 += ".fam";
	int N1 = getLineNum(famfile1);
	string bimfile1 = stringname1;
	bimfile1 += ".bim";
	long int P1 = getLineNum(bimfile1);

	IntegerVector A11(P1), A12(P1);
	CharacterVector rsname1(P1);
	IntegerVector chr1(P1), bp1(P1);
	NumericVector morgan1(P1);

	//read from SNPinfo file (pass as pointer)
	cout << "## Start loading SNP info:" << endl;
	clock_t t1 = clock();
	ReadSNPinfo(bimfile1, A11, A12, rsname1, chr1, bp1, morgan1, P1);
	cout << "Finish loading SNP info in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;

	// load IID in file 1 (fam file 1)
	cout << "## Start loading fam files 1 " << endl;
	IntegerVector tmp1(N1);
	CharacterVector FID_1(N1), IID_1(N1);
	NumericVector tmp2(N1);
	ReadPlinkFamFile(famfile1, FID_1, IID_1, tmp1, tmp2, N1);
	cout << "Finish loading fam files 1 in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;

	ObjXY obj1_XY = ReadDataFromFile(stringname1);
	Mat<unsigned>* X1_point = &obj1_XY.X;


	//-------------------------------------------------------------------------------------------//
	cout << "begin to read GWAS data" << endl;

	int P2 = getLineNum(stringname2);
	List GWAS = read_GWAS(stringname2, P2);

	IntegerVector A21(P2), A22(P2);
	CharacterVector rsname2(P2);
	IntegerVector chr2(P2), bp2(P2);

	A21 = GWAS["A1"];
	A22 = GWAS["A2"];
	rsname2 = GWAS["GWAS_snps"];
	chr2 = GWAS["GWAS_chr"];
	bp2 = GWAS["BP"];
	mat GWAS_sum = GWAS["GWAS_sum"];


	//-------------------------------------------------------------------------------------------//
	cout << "begin to read reference panel" << endl;

	string famfile3 = stringname3;
	famfile3 += ".fam";
	int N3 = getLineNum(famfile3);
	string bimfile3 = stringname3;
	bimfile3 += ".bim";
	long int P3 = getLineNum(bimfile3);

	IntegerVector A31(P3), A32(P3);
	CharacterVector rsname3(P3);
	IntegerVector chr3(P3), bp3(P3);
	NumericVector morgan(P3);

	//read from SNPinfo file (pass as pointer)
	cout << endl;
	cout << "Start loading SNP info:" << endl;
	clock_t t3 = clock();
	ReadSNPinfo(bimfile3, A31, A32, rsname3, chr3, bp3, morgan, P3);
	cout << "Finish loading SNP info in " << (clock() - t3)*1.0 / CLOCKS_PER_SEC << " sec." << endl;

	ObjXY obj3_XY = ReadDataFromFile(stringname3);
	Mat<unsigned>* X3_point = &obj3_XY.X;


	//-------------------------------------------------------------------------------------------//
	cout << "## Start loading expression file. " << endl;
	int Ngene = getLineNum(stringname4) - 1; // header is first line;
	char delimiter = '\t';
	List tmp = getColNum_Header2(stringname4, delimiter);

	int Nindiv = tmp["columns"];
	CharacterVector header = tmp["header"];

	// IID from header of expr file
	uvec idxtmp = linspace<uvec>(6, Nindiv - 1, Nindiv - 6);
	CharacterVector IID_expr = charv_subset2_(header, idxtmp);

	// match file 1 with expression file
	CharacterVector indiv_inter = intersect(IID_1, IID_expr);
	IntegerVector idxint = match(IID_1, indiv_inter); // put the IID in the order of file 1
	CharacterVector indiv_4use = IID_1[Rcpp::is_na(idxint) == false];
	IntegerVector idxin1 = match(indiv_4use, IID_1) - 1;  //index for individual in file 1
	IntegerVector idxin2 = match(indiv_4use, IID_expr) - 1;  //index for individual in expr file

	uvec idxin1_ = as<uvec>(idxin1);
	uvec idxin2_ = as<uvec>(idxin2);

	//read expr file
	List expr_tmp = getExprFile2(stringname4, delimiter, Nindiv, Ngene + 1, 6);
	CharacterVector genetype1 = expr_tmp["genetype1"];
	CharacterVector genetype2 = expr_tmp["genetype2"];
	CharacterVector targetID = expr_tmp["targetID"];
	mat expr = expr_tmp["expr"];
	vec lower = expr_tmp["lower"];
	vec upper = expr_tmp["upper"];
	vec chr_expr = expr_tmp["chr_expr"];

	// subset overlapped individuals
	mat* exprTmp;
	mat expr_used;
	exprTmp = &expr;
	expr_used = exprTmp->cols(as<uvec>(idxin2));

	//-------------------------------------------------------------------------------------------//
	cout << "matching SNPs and correcting the direction of minor allele" << endl;

	CharacterVector rs_12 = intersect(rsname1, rsname2);
	CharacterVector rs_123 = intersect(rs_12, rsname3);

	IntegerVector idxin = match(rsname1, rs_123);
	CharacterVector rsname_4use = rsname1[Rcpp::is_na(idxin) == false]; // rsname_4use in the order of eQTL.
	IntegerVector bp1_4use = bp1[Rcpp::is_na(idxin) == false]; // bp1 in the order of eQTL
	IntegerVector chr1_4use = chr1[Rcpp::is_na(idxin) == false]; // chr1 in the order of eQTL

	IntegerVector idx_in1 = match(rsname_4use, rsname1);  //index for SNPs in eQTL
	IntegerVector idx_in2 = match(rsname_4use, rsname2);  //index for SNPs in GWAS
	IntegerVector idx_in3 = match(rsname_4use, rsname3);  //index for SNPs in reference panel
	cout << "Size of matched SNPs: " << rsname_4use.size() << endl;

	//----------------------------------------------------------------------------------------------//
	//compare direction
	uvec temp;

	uvec idx1 = as<uvec>(idx_in1) -1;
	temp = as<uvec>(A11);
	uvec A11_ = temp.elem(idx1);
	temp = as<uvec>(A12);
	uvec A12_ = temp.elem(idx1);

	uvec idx2 = as<uvec>(idx_in2) -1;
	temp = as<uvec>(A21);
	uvec A21_ = temp.elem(idx2);
	temp = as<uvec>(A22);
	uvec A22_ = temp.elem(idx2);

	uvec idx3 = as<uvec>(idx_in3) -1;
	temp = as<uvec>(A31);
	uvec A31_ = temp.elem(idx3);
	temp = as<uvec>(A32);
	uvec A32_ = temp.elem(idx3);

	//----------------------------------------------------------------------------------------------//
	//ascii: (A:65; C:67; G:71; T:84) (a:97;c:99,g:103;t:116)
	//replace lower letter to upper letter in A11_, A12_, A21_ and A22_;

	uvec idx;

	idx = find(A11_ > 85);
	A11_.elem(idx) = A11_.elem(idx) - 32;
	idx = find(A12_ > 85);
	A12_.elem(idx) = A12_.elem(idx) - 32;
	idx = find(A11_ > 85);
	A21_.elem(idx) = A21_.elem(idx) - 32;
	idx = find(A12_ > 85);
	A22_.elem(idx) = A22_.elem(idx) - 32;
	idx = find(A31_ > 85);
	A31_.elem(idx) = A31_.elem(idx) - 32;
	idx = find(A32_ > 85);
	A32_.elem(idx) = A32_.elem(idx) - 32;


	//A11_: replace T with A,replace G with C
	T2AG2C(A11_);

	//A12_: replace T with A,replace G with C
	T2AG2C(A12_);

	//A31_: replace T with A, replace G with C
	T2AG2C(A31_);

	//A32_: replace T with A,replace G with C
	T2AG2C(A32_);

	//A21_: replace T with A,replace G with C
	T2AG2C(A21_);

	//A22_: replace T with A,replace G with C
	T2AG2C(A22_);

	//remove index having the same allele
	uvec idxtmp1, idxtmp2;
	idxtmp1 = find((A31_ + A32_) == (A11_ + A12_));
	idxtmp2 = find((A31_ + A32_) == (A21_ + A22_));

	uvec idx4 = intersect(idxtmp1, idxtmp2);

	uvec A11_r = A11_.elem(idx4), A12_r = A12_.elem(idx4);
	uvec A21_r = A21_.elem(idx4), A22_r = A22_.elem(idx4);
	uvec A31_r = A31_.elem(idx4), A32_r = A32_.elem(idx4);

	IntegerVector idx_tmp = wrap(idx4);
	CharacterVector rsname_4use_r = rsname_4use[idx_tmp];
	IntegerVector bp_4use_r = bp1_4use[idx_tmp];
	IntegerVector chr_4use_r = chr1_4use[idx_tmp];

	cout << "Size of A31_r is " << A31_r.n_elem << endl;
	cout << "Size of A11_r is " << A11_r.n_elem << endl;
	cout << "Size of A21_r is " << A21_r.n_elem << endl;

	IntegerVector idxinFile1_ = match(rsname_4use_r, rsname1) - 1;  //index for SNPs in eQTL
	IntegerVector idxinFile2_ = match(rsname_4use_r, rsname2) - 1;  //index for SNPs in GWAS
	IntegerVector idxinFile3_ = match(rsname_4use_r, rsname3) - 1; //index for SNPs in reference panel

	uvec idxinFile1 = as<uvec>(idxinFile1_);
	uvec idxinFile2 = as<uvec>(idxinFile2_);
	uvec idxinFile3 = as<uvec>(idxinFile3_);

	uvec idx21 = find(A11_r != A21_r);
	cout << "Size of A11_r != A21_r is " << idx21.n_elem << endl;
	uvec temp1 = idx2(idx4);
	GWAS_sum(temp1(idx21), zeros<uvec>(1)) = -GWAS_sum(temp1(idx21), zeros<uvec>(1));

	uvec idx22 = find(A11_r != A31_r);
	cout << "Size of A11_r != A31_r is " << idx22.n_elem << endl;
	uvec temp2 = idx3(idx4);
	X3_point->cols(temp2(idx22)) = 2 - X3_point->cols(temp2(idx22));

	//----------------------------------------------------------------------------------------------//

	mat expr_info = zeros<mat>(lower.n_elem, 3);
	expr_info.col(0) = chr_expr;
	expr_info.col(1) = lower;
	expr_info.col(2) = upper;

	vec bp_4use_r_vec = as<vec>(bp_4use_r);
	vec chr_4use_r_vec = as<vec>(chr_4use_r);

	mat snp_info = zeros<mat>(chr_4use_r.length(), 2);
	snp_info.col(0) = bp_4use_r_vec;
	snp_info.col(1) = chr_4use_r_vec;

	//---------------------------------------------------------------//
	// centering
	int size2 = X1_point->n_cols;
	vec meanX(size2);
	mat X1 = zeros(N1, size2);

	for (int i = 0; i < size2; i++){
		meanX[i] = sum(X1_point->col(i))*1.0 / N1;
		vec v_i = conv_to<vec>::from(X1_point->col(i));
		v_i = v_i - meanX[i];
		X1.col(i) = v_i;
	}

	int size3 = X3_point->n_cols;
	vec meanX3(size3);
	mat X3 = zeros(N3, size3);

	for (int i = 0; i < size3; i++){
		meanX3[i] = sum(X3_point->col(i))*1.0 / N3;
		vec v_i = conv_to<vec>::from(X3_point->col(i));
		v_i = v_i - meanX3[i];
		X3.col(i) = v_i;
	}

	//----------------------------------------------------------------//
	// testing
	t1 = clock();
	cout << "### Precompute the active genes ... ";
	uvec idx_active_gene = zeros<uvec>(0);
	uvec g_tmp(1);
	for (uword g = 0; g < Ngene; g++){

		idx = find(bp_4use_r_vec < upper(g) + bw && bp_4use_r_vec > lower(g) - bw
			&& chr_4use_r_vec == chr_expr(g));
		if (idx.is_empty() == false){
			if (idx.n_elem > 1){
				g_tmp(0) = g;
				idx_active_gene = join_cols(idx_active_gene, g_tmp);
			}
		}
	}
	cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;


	// combine gene info into mat
	uword Ngene_active = idx_active_gene.n_elem;
	mat gene_info = zeros<mat>(Ngene_active, 3);
	gene_info = expr_info.rows(idx_active_gene);
	mat expr_active = expr_used.rows(idx_active_gene);

	CharacterVector gene_type1 = genetype1[as<Rcpp::IntegerVector>(wrap(idx_active_gene))];
	CharacterVector gene_type2 = genetype2[as<Rcpp::IntegerVector>(wrap(idx_active_gene))];

	vec stat = zeros(Ngene_active, 1);
	vec alphag = zeros(Ngene_active, 1);
	uvec NumSNP = zeros<uvec>(Ngene_active, 1);

	cout << "### Start running CoMM_S2 on active genes ... " << endl;;
	t1 = clock();

	mat out_param0 = -99 * ones<mat>(Ngene_active, 7);

	for (int g = 0; g < Ngene_active; g++){
		uvec idx = find(snp_info.col(0) < gene_info(g, 2) + bw && snp_info.col(0) > gene_info(g, 1) - bw
			&& snp_info.col(1) == gene_info(g, 0));
		if (idx.n_elem == 0){
			cout << "Error: no matching SNPs for " << g << "-th gene ... " << endl;
		}

		// extract the SNPs within the region as X1 and X3
		mat X1tmp = X1.cols(idxinFile1(idx));
		X1tmp = X1tmp.rows(idxin1_);
		mat X3tmp = X3.cols(idxinFile3(idx));

		uword idx1 = g;
		vec y = trans(expr_active.row(idx1));

		// remove the missing individuals from eQTL data
		uvec idx2 = find_finite(y);
		X1tmp = X1tmp.rows(idx2);
		y = y(idx2);

		// calculate cellular heritability
		rowvec sdX1tmp = stddev(X1tmp, 0, 0); // see manual
		mat X1tmp2 = X1tmp / repmat(sdX1tmp, X1tmp.n_rows, 1);

		// attention 
		uvec idx_finite = find_finite(X1tmp2.row(0));
		X1tmp2 = X1tmp2.cols(idx_finite);
		X1tmp = X1tmp.cols(idx_finite);
		X3tmp = X3tmp.cols(idx_finite);

		int Mg = X3tmp.n_cols;

		cout << "Info: the number of snps for " << g << "-th gene is " << Mg << endl;

		uvec idx_col1 = zeros<uvec>(1);
		uvec idx_col2 = ones<uvec>(1);
		vec hatmu2 = GWAS_sum(idxinFile2(idx), idx_col1);
		vec hats2 = GWAS_sum(idxinFile2(idx), idx_col2);

		mat R1 = eye(Mg, Mg);

		vec diagxtx3 = sum(X3tmp%X3tmp).t();
		for (int i = 0; i < Mg; i++){
			for (int j = 0; j < Mg; j++){
				if (i != j){
					R1(i, j) = sum(X3tmp.col(i) % X3tmp.col(j)) / sqrt(diagxtx3(i)*diagxtx3(j));
				}
			}
		}
		mat R;

		// positive definite sparse correlation matrix
		R = lam*R1 + (1 - lam)*eye(Mg, Mg);


		// set the parameter
		Options_SCLMM* lp_opt = new Options_SCLMM(100000, 0, 10, 1e-5, 0);
		Options_SCLMM* lp_opt1 = new Options_SCLMM(100000, 0, 10, 1e-5, 1);

		mat w1tmp = ones(X1tmp.n_rows, 1);

		double sigma2y, sigma2beta, loglik;
		vec beta0 = zeros<vec>(w1tmp.n_cols);
		int iterLmm, maxIter = 100000;
		mat Sigb = zeros<mat>(X1tmp.n_cols, X1tmp.n_cols);
		vec mub = zeros<vec>(X1tmp.n_cols);

		lmm_pxem_ptr(y, w1tmp, X1tmp2, maxIter, sigma2y, sigma2beta, beta0, loglik, iterLmm, Sigb, mub);

		out_param0(g, 0) = sigma2beta*X1tmp.n_cols;
		out_param0(g, 1) = sigma2y;
		out_param0(g, 6) = 1 / (1 + sigma2y / (sigma2beta*X1tmp.n_cols));

		ObjSCLMM objHa = rcpparma_SCLMM_IS(X1tmp, y, w1tmp, hatmu2, hats2, R, lp_opt, 1);
		ObjSCLMM objH0 = rcpparma_SCLMM_IS(X1tmp, y, w1tmp, hatmu2, hats2, R, lp_opt1, 1);

		alphag(g) = objHa.alphag;
		stat(g) = 2 * (objHa.LRLB - objH0.LRLB);
		NumSNP(g) = Mg;

		out_param0(g, 3) = 2 * (objHa.LRLB - objH0.LRLB);
		out_param0(g, 4) = objHa.alphag;
		out_param0(g, 5) = Mg;

		// reset
		X1tmp.reset();
		X3tmp.reset();
		hatmu2.reset();
		hats2.reset();


		if ((g + 1) % 100 == 0 && (g + 1) != 0){
			cout << g + 1 << "-th Gene starts working ..." << endl;
			cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
		}

	}

	mat out_param = -99 * ones<mat>(Ngene, 7);
	out_param.rows(idx_active_gene) = out_param0;


	List output = List::create(
		Rcpp::Named("expr") = expr,  // expression
		Rcpp::Named("expr_active") = expr_active, //expression for active
		Rcpp::Named("rsname_4use_r") = rsname_4use_r,
		Rcpp::Named("snp_info_4use_r") = snp_info,
		Rcpp::Named("targetID") = targetID,
		Rcpp::Named("genetype1") = genetype1,
		Rcpp::Named("genetype2") = genetype2, // all gene name
		Rcpp::Named("gene_type1") = gene_type1,
		Rcpp::Named("gene_type2") = gene_type2, // active gene name
		Rcpp::Named("expr_info") = expr_info, // gene information for all genes: chr lower upper
		Rcpp::Named("gene_info") = gene_info, // gene information for active genes: chr lower upper
		Rcpp::Named("stat") = stat,
		Rcpp::Named("alphag") = alphag,
		Rcpp::Named("NumSNP") = NumSNP,
		Rcpp::Named("out_param") = out_param,
		Rcpp::Named("out_param0") = out_param0);
	return output;

}





