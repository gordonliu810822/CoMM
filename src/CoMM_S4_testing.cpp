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
#include "CoMM_S4.hpp"
#include <boost/algorithm/string.hpp>
#include <R.h>
#include "readPlink2.hpp"
#include "StandardizeData.hpp"
#include "plinkfun.hpp"
#include "mt_paral_testing_job3.hpp"
#include "CoMM_S2_testing.hpp"

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





// [[Rcpp::export]]
List read_eqtls(std::string filename, int P){
  std::ifstream ifs(filename.c_str());
  std::string line;
  vector<string> fields;
  CharacterVector snp(P);
  CharacterVector gene_all(P);
  mat eqtls_sum = zeros(P, 2); // transpcriptome summary data
  int i = 0;
  CharacterVector gene_(1);
  while (std::getline(ifs, line)) // read one line from ifs
  {
    std::istringstream iss(line); // access line as a stream
    boost::split(fields, line, boost::is_any_of(" \t *"));
    
    snp[i] = fields[0];
    gene_all[i] = fields[1];
    eqtls_sum(i, 0) = atof(fields[2].c_str());
    eqtls_sum(i, 1) = atof(fields[3].c_str());
    i++;
  }
  ifs.close();
     
  // Record gene information   
  CharacterVector gene = unique(gene_all);   
  IntegerVector gene_first_occur_idx = match(gene, gene_all);
  int Ngene = gene.length();
  umat gene_index = zeros<umat>(Ngene, 2);
  IntegerVector tmp = sort_unique(gene_first_occur_idx);
  IntegerVector idx =  match(tmp, gene_first_occur_idx);
  CharacterVector gene2(Ngene);
  for (int i = 0; i < Ngene; i++){
    if (i < Ngene - 1){
      gene_index(i,0) = tmp[i];
      gene_index(i,1) = tmp[i+1] - 1;
      gene2[i] = gene[idx[i]-1];
    }
    if (i == Ngene - 1){
      gene_index(i,0) = tmp[i];
      gene_index(i,1) = P;
      gene2[i] = gene[idx[i]-1];
      }
  }
  //cout << gene << endl;
  //cout << idx << endl;
   
  List output = List::create(
    Rcpp::Named("eqtls_snps") = snp,
    Rcpp::Named("eqtls_gene") = gene2,
    Rcpp::Named("eqtls_sum") = eqtls_sum,
    Rcpp::Named("gene_index") = gene_index);
  
  return output;

}


// [[Rcpp::export]]
List CoMM_S4_testing(std::string stringname1, std::string stringname2, std::string stringname3, std::string stringname4, std::string stringname5, bool px = 1, double lam = 0.95){
	// Summary data for transpcriptome data: stringname1; 
	// Summary data for GWAS data; stringname2;
	// Individual data for eQTL reference panel; stringname3;
	// Individual data for GWAS reference panel; stringname5;

	int P1 = getLineNum(stringname1);
	int P2 = getLineNum(stringname2);

	string famfile3 = stringname3;
	famfile3 += ".fam";
	int N3 = getLineNum(famfile3);
	string bimfile3 = stringname3;
	bimfile3 += ".bim";
	long int P3 = getLineNum(bimfile3);

	string famfile5 = stringname5;
	famfile5 += ".fam";
	int N5 = getLineNum(famfile5);
	string bimfile5 = stringname5;
	bimfile5 += ".bim";
	long int P5 = getLineNum(bimfile5);

	Rcout << "begin to read reference panel" << endl;
	ObjXY obj_XY = ReadDataFromFile(stringname3);
	arma::Mat<unsigned>* X3_point = &obj_XY.X;
	ObjXY obj2_XY = ReadDataFromFile(stringname5);
	arma::Mat<unsigned>* X5_point = &obj2_XY.X;


	string bimfile4 = stringname4;
	bimfile4 += ".bim";
	long int P4 = getLineNum(bimfile4);

	List eqtls = read_eqtls(stringname1, P1);
	List GWAS = read_GWAS(stringname2, P2);

	CharacterVector eqtls_snp = eqtls["eqtls_snps"];
	CharacterVector eqtls_gene = eqtls["eqtls_gene"];
	mat eqtls_sum = eqtls["eqtls_sum"];
	umat gene_index = eqtls["gene_index"];

	CharacterVector GWAS_snp = GWAS["GWAS_snps"];
	mat GWAS_sum = GWAS["GWAS_sum"];


	////////////////////  eQTL SNP information  //////////////////////////////////////////////
	IntegerVector A11(P4), A12(P4);
	CharacterVector rsname1(P4);
	IntegerVector chr1(P4), bp1(P4);
	NumericVector morgan1(P4);

	//read from SNPinfo file (pass as pointer)
	cout << "Start loading SNP info:" << endl;
	clock_t t1 = clock();
	ReadSNPinfo(bimfile4, A11, A12, rsname1, chr1, bp1, morgan1, P4);
	cout << "Finish loading SNP info in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;


	////////////////////  GWAS SNP information  //////////////////////////////////////////////
	IntegerVector A21(P2), A22(P2);
	CharacterVector rsname2(P2);
	IntegerVector chr2(P2), bp2(P2);

	A21 = GWAS["A1"];
	A22 = GWAS["A2"];
	rsname2 = GWAS["GWAS_snps"];
	chr2 = GWAS["GWAS_chr"];
	bp2 = GWAS["BP"];


	//////////////////// First reference panel SNP information   //////////////////////////////
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


	//////////////////// Second reference panel SNP information   //////////////////////////////
	IntegerVector A51(P5), A52(P5);
	CharacterVector rsname5(P5);
	IntegerVector chr5(P5), bp5(P5);
	NumericVector morgan5(P5);

	//read from SNPinfo file (pass as pointer)
	cout << endl;
	cout << "Start loading SNP info:" << endl;
	clock_t t5 = clock();
	ReadSNPinfo(bimfile5, A51, A52, rsname5, chr5, bp5, morgan5, P5);
	cout << "Finish loading SNP info in " << (clock() - t5)*1.0 / CLOCKS_PER_SEC << " sec." << endl;


	///////////////////  mathcing SNP  ////////////////////////////////////////////////////////////
	CharacterVector rs_12 = intersect(rsname1, rsname2);
	CharacterVector rs_123 = intersect(rs_12, rsname3);
	CharacterVector rs_1235 = intersect(rs_123, rsname5);
	//cout << "rsname1 " << rsname1.length() << endl;
	//cout << "rsname2 " << rsname2.length() << endl;
	//cout << "rsname3 " << rsname3.length() << endl;
	//cout << "rsname5 " << rsname5.length() << endl;
	//cout << "rs_12 " << rs_12.length() << endl;
	//cout << "rs_123 " << rs_123.length() << endl;
	//cout << "rs_1235 " << rs_1235.length() << endl;


	//////////////////  rsname_4use ///////////////////////////////////////////////////////////////
	IntegerVector idxin = match(rsname3, rs_1235);
	CharacterVector rsname_4use = rsname3[Rcpp::is_na(idxin) == false]; // rsname_4use in the order of panel.
	IntegerVector bp3_4use = bp3[Rcpp::is_na(idxin) == false]; // bp3 in the order of panel
	IntegerVector chr3_4use = chr3[Rcpp::is_na(idxin) == false]; // chr3 in the order of panel

	//----------------------------------------------------------------------------------------------//
	printf("Start matching SNPs with summary stat for eqtls. \n");
	IntegerVector idxin1 = match(rsname_4use, rsname1);  //index for SNPs in eqtls

	//----------------------------------------------------------------------------------------------//
	printf("Start matching SNPs with summary stat for GWAS. \n");
	IntegerVector idxin2 = match(rsname_4use, rsname2);  //index for SNPs in GWAS

	//----------------------------------------------------------------------------------------------//
	printf("Start matching SNPs with first reference panel. \n");
	IntegerVector idxin3 = match(rsname_4use, rsname3); //index for SNPs in panel SNPs

	//----------------------------------------------------------------------------------------------//
	printf("Start matching SNPs with second reference panel. \n");
	IntegerVector idxin5 = match(rsname_4use, rsname5); //index for SNPs in panel SNPs

	cout << "Size of matched SNPs: " << rsname_4use.size() << endl; //compare keepIndx in R: Height_1000Q_match.R


	//////////////////   correct direction of minor allele  ////////////////////////////////////////////
	uvec tmp;

	uvec idx1 = as<uvec>(idxin1) -1;
	tmp = as<uvec>(A11);
	uvec A11_ = tmp.elem(idx1);
	tmp = as<uvec>(A12);
	uvec A12_ = tmp.elem(idx1);

	uvec idx2 = as<uvec>(idxin2) -1;
	tmp = as<uvec>(A21);
	uvec A21_ = tmp.elem(idx2);
	tmp = as<uvec>(A22);
	uvec A22_ = tmp.elem(idx2);

	uvec idx3 = as<uvec>(idxin3) -1;
	tmp = as<uvec>(A31);
	uvec A31_ = tmp.elem(idx3);
	tmp = as<uvec>(A32);
	uvec A32_ = tmp.elem(idx3);

	uvec idx5 = as<uvec>(idxin5) -1;
	tmp = as<uvec>(A51);
	uvec A51_ = tmp.elem(idx5);
	tmp = as<uvec>(A52);
	uvec A52_ = tmp.elem(idx5);

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
	idx = find(A51_ > 85);
	A51_.elem(idx) = A51_.elem(idx) - 32;
	idx = find(A52_ > 85);
	A52_.elem(idx) = A52_.elem(idx) - 32;


	//A11_: replace T with A,replace G with C
	T2AG2C(A11_);

	//A12_: replace T with A,replace G with C
	T2AG2C(A12_);

	//A21_: replace T with A,replace G with C
	T2AG2C(A21_);

	//A22_: replace T with A,replace G with C
	T2AG2C(A22_);

	//A31_: replace T with A, replace G with C
	T2AG2C(A31_);

	//A32_: replace T with A,replace G with C
	T2AG2C(A32_);

	//A51_: replace T with A,replace G with C
	T2AG2C(A51_);

	//A52_: replace T with A,replace G with C
	T2AG2C(A52_);

	//remaining index having the same allele
	uvec idxtmp1, idxtmp2, idxtmp5;
	idxtmp1 = find((A31_ + A32_) == (A11_ + A12_));
	idxtmp2 = find((A31_ + A32_) == (A21_ + A22_));
	idxtmp5 = find((A31_ + A32_) == (A51_ + A52_));
	//cout << A11_.n_elem << "  " << A21_.n_elem << "  " << A31_.n_elem << "  " << A51_.n_elem << endl;
	//cout << idxtmp1.n_elem << "   " << idxtmp2.n_elem << "   " << idxtmp5.n_elem << endl;
	uvec idx4_tmp = intersect(idxtmp1, idxtmp2);
	uvec idx4 = intersect(idx4_tmp, idxtmp5);

	uvec A11_r = A11_.elem(idx4), A12_r = A12_.elem(idx4);
	uvec A21_r = A21_.elem(idx4), A22_r = A22_.elem(idx4);
	uvec A31_r = A31_.elem(idx4), A32_r = A32_.elem(idx4);
	uvec A51_r = A51_.elem(idx4), A52_r = A52_.elem(idx4);

	IntegerVector idxtmp = wrap(idx4);
	CharacterVector rsname_4use_r = rsname_4use[idxtmp];
	IntegerVector bp_4use_r = bp3_4use[idxtmp];
	IntegerVector chr_4use_r = chr3_4use[idxtmp];

	cout << "Size of A11_r is " << A11_r.n_elem << endl;
	cout << "Size of A21_r is " << A21_r.n_elem << endl;
	cout << "Size of A31_r is " << A31_r.n_elem << endl;
	cout << "Size of A51_r is " << A51_r.n_elem << endl;

	uvec idx21 = find(A11_r != A21_r);
	cout << "Size of SNPs with opposite allele direction for GWAS compare with eQTL " << idx21.n_elem << endl;
	uvec tmp1 = idx2(idx4);
	GWAS_sum(tmp1(idx21), zeros<uvec>(1)) = -GWAS_sum(tmp1(idx21), zeros<uvec>(1));

	uvec idx22 = find(A11_r != A31_r);
	cout << "Size of SNPs with opposite allele direction for First reference panel compare with eQTL " << idx22.n_elem << endl;
	uvec tmp2 = idx3(idx4);
	X3_point->cols(tmp2(idx22)) = 2 - X3_point->cols(tmp2(idx22));

	uvec idx25 = find(A11_r != A51_r);
	cout << "Size of SNPs with opposite allele direction for Second reference panel compare with eQTL " << idx25.n_elem << endl;
	uvec tmp5 = idx5(idx4);
	X5_point->cols(tmp5(idx25)) = 2 - X5_point->cols(tmp5(idx25));


	///////////////////  begin to fit the model //////////////////////////
	mat X = StandardX(X3_point);
	mat X2 = StandardX(X5_point);

	// set the parameter
	Options_CoMM_S4* lp_opt = new Options_CoMM_S4(10000, 0, 10, 1e-5, 0);
	Options_CoMM_S4* lp_opt1 = new Options_CoMM_S4(10000, 0, 10, 1e-5, 1);

	// number of gene
	int Ngene = gene_index.n_rows;
	vec stat = zeros(Ngene, 1);
	vec alphag = zeros(Ngene, 1);
	uvec NumSNP = zeros<uvec>(Ngene, 1);


	uvec idx_all = zeros<uvec>(0);
	uvec idx_active_gene = zeros<uvec>(0);
	uvec g_tmp(1);


	// testing 
	cout << "### Start running CoMM_2S2 on active genes ... " << endl;
	t1 = clock();
	mat out_param = -99 * ones<mat>(Ngene, 3);

	for (int g = 0; g < Ngene; g++){
		if (g % 100 == 0 && g != 0){
			cout << g << "-th Gene starts working ...";
			cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
		}

		int Nsnp_gene_eQTL = gene_index(g, 1) - gene_index(g, 0) + 1;
		IntegerVector idx_chr = wrap(linspace<uvec>(gene_index(g, 0), gene_index(g, 1), Nsnp_gene_eQTL) - 1);
		CharacterVector snp_gene_eqtls = eqtls_snp[idx_chr];

		IntegerVector idx_4use_r = match(snp_gene_eqtls, rsname_4use_r);
		CharacterVector snp_gene = snp_gene_eqtls[Rcpp::is_na(idx_4use_r) == false]; //snp in the order of snp_each_gene

		int Nsnp_gene = snp_gene.length();

		if (Nsnp_gene > 1){
			g_tmp(0) = g;
			idx_active_gene = join_cols(idx_active_gene, g_tmp);

			IntegerVector index1 = match(snp_gene, snp_gene_eqtls);
			uvec idx31 = as<uvec>(index1) -1;
			IntegerVector index2 = match(snp_gene, rsname2);
			uvec idx32 = as<uvec>(index2) -1;
			IntegerVector index3 = match(snp_gene, rsname3);
			uvec idx33 = as<uvec>(index3) -1;
			IntegerVector index5 = match(snp_gene, rsname5);
			uvec idx35 = as<uvec>(index5) -1;

			uvec idx_col1 = zeros<uvec>(1);
			uvec idx_col2 = ones<uvec>(1);

			uvec idx41 = linspace<uvec>(gene_index(g, 0), gene_index(g, 1), Nsnp_gene_eQTL) - 1;
			uvec idx51 = idx41(idx31);

			vec hatmu1 = eqtls_sum(idx51, idx_col1);
			vec hats1 = eqtls_sum(idx51, idx_col2);
			vec hatmu2 = GWAS_sum(idx32, idx_col1);
			vec hats2 = GWAS_sum(idx32, idx_col2);
			vec zscore1 = hatmu1 / hats1;
			vec zscore2 = hatmu2 / hats2;

			cout << "Info: the number of snps for " << g << "-th gene is " << Nsnp_gene << endl;

			//First reference panel
			mat X_sub = X.cols(idx33);
			mat R_emp = eye(Nsnp_gene, Nsnp_gene);
			for (int i = 0; i < Nsnp_gene; i++){
				for (int j = 0; j < Nsnp_gene; j++){
					if (i != j){
						R_emp(i, j) = sum(X_sub.col(i) % X_sub.col(j)) / (norm(X_sub.col(i))*norm(X_sub.col(j)));
					}
				}
			}
			mat R = lam*R_emp + (1 - lam)*eye(R_emp.n_cols, R_emp.n_cols);

			//Second reference panel
			mat X2_sub = X2.cols(idx35);
			mat R2_emp = eye(Nsnp_gene, Nsnp_gene);
			for (int i = 0; i < Nsnp_gene; i++){
				for (int j = 0; j < Nsnp_gene; j++){
					if (i != j){
						R2_emp(i, j) = sum(X2_sub.col(i) % X2_sub.col(j)) / (norm(X2_sub.col(i))*norm(X2_sub.col(j)));
					}
				}
			}
			mat R2 = lam*R2_emp + (1 - lam)*eye(R2_emp.n_cols, R2_emp.n_cols);


			ObjCoMM_S4 objHa = rcpparma_CoMM_S4_VB(zscore1, zscore2, R, R2, lp_opt, px);
			ObjCoMM_S4 objH0 = rcpparma_CoMM_S4_VB(zscore1, zscore2, R, R2, lp_opt1, px);

			alphag(g) = objHa.alphag;
			stat(g) = 2 * (objHa.LRLB - objH0.LRLB);
			NumSNP(g) = index2.length();

			out_param(g, 0) = 2 * (objHa.LRLB - objH0.LRLB);
			out_param(g, 1) = objHa.alphag;
			out_param(g, 2) = index2.length();

			// reset
			X_sub.reset();
			X2_sub.reset();
			hatmu1.reset();
			hats1.reset();
			hatmu2.reset();
			hats2.reset();

		}

	}

	mat out_param0 = out_param.rows(idx_active_gene);
	CharacterVector gene_type1 = eqtls_gene[as<Rcpp::IntegerVector>(wrap(idx_active_gene))];


	List output = List::create(
		Rcpp::Named("genetype1") = eqtls_gene,
		Rcpp::Named("gene_type1") = gene_type1,
		Rcpp::Named("out_param") = out_param,
		Rcpp::Named("out_param0") = out_param0);
	return output;

}



// [[Rcpp::export]]
List CoMM_S4_testing_mt(std::string stringname1, std::string stringname2, std::string stringname3, std::string stringname4, std::string stringname5, bool px = 1, double lam = 0.95, int coreNum = 24){
	// Summary data for transpcriptome data: stringname1; 
	// Summary data for GWAS data; stringname2;
	// Individual data for eQTL reference panel; stringname3;
	// Individual data for GWAS reference panel; stringname5;

	int P1 = getLineNum(stringname1);
	int P2 = getLineNum(stringname2);

	string famfile3 = stringname3;
	famfile3 += ".fam";
	int N3 = getLineNum(famfile3);
	string bimfile3 = stringname3;
	bimfile3 += ".bim";
	long int P3 = getLineNum(bimfile3);

	string famfile5 = stringname5;
	famfile5 += ".fam";
	int N5 = getLineNum(famfile5);
	string bimfile5 = stringname5;
	bimfile5 += ".bim";
	long int P5 = getLineNum(bimfile5);

	Rcout << "begin to read reference panel" << endl;
	ObjXY obj_XY = ReadDataFromFile(stringname3);
	arma::Mat<unsigned>* X3_point = &obj_XY.X;
	ObjXY obj2_XY = ReadDataFromFile(stringname5);
	arma::Mat<unsigned>* X5_point = &obj2_XY.X;


	string bimfile4 = stringname4;
	bimfile4 += ".bim";
	long int P4 = getLineNum(bimfile4);

	List eqtls = read_eqtls(stringname1, P1);
	List GWAS = read_GWAS(stringname2, P2);

	CharacterVector eqtls_snp = eqtls["eqtls_snps"];
	CharacterVector eqtls_gene = eqtls["eqtls_gene"];
	mat eqtls_sum = eqtls["eqtls_sum"];
	umat gene_index = eqtls["gene_index"];

	CharacterVector GWAS_snp = GWAS["GWAS_snps"];
	mat GWAS_sum = GWAS["GWAS_sum"];


	////////////////////  eQTL SNP information  //////////////////////////////////////////////
	IntegerVector A11(P4), A12(P4);
	CharacterVector rsname1(P4);
	IntegerVector chr1(P4), bp1(P4);
	NumericVector morgan1(P4);

	//read from SNPinfo file (pass as pointer)
	cout << "Start loading SNP info:" << endl;
	clock_t t1 = clock();
	ReadSNPinfo(bimfile4, A11, A12, rsname1, chr1, bp1, morgan1, P4);
	cout << "Finish loading SNP info in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;


	////////////////////  GWAS SNP information  //////////////////////////////////////////////
	IntegerVector A21(P2), A22(P2);
	CharacterVector rsname2(P2);
	IntegerVector chr2(P2), bp2(P2);

	A21 = GWAS["A1"];
	A22 = GWAS["A2"];
	rsname2 = GWAS["GWAS_snps"];
	chr2 = GWAS["GWAS_chr"];
	bp2 = GWAS["BP"];


	//////////////////// First reference panel SNP information   //////////////////////////////
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


	//////////////////// Second reference panel SNP information   //////////////////////////////
	IntegerVector A51(P5), A52(P5);
	CharacterVector rsname5(P5);
	IntegerVector chr5(P5), bp5(P5);
	NumericVector morgan5(P5);

	//read from SNPinfo file (pass as pointer)
	cout << endl;
	cout << "Start loading SNP info:" << endl;
	clock_t t5 = clock();
	ReadSNPinfo(bimfile5, A51, A52, rsname5, chr5, bp5, morgan5, P5);
	cout << "Finish loading SNP info in " << (clock() - t5)*1.0 / CLOCKS_PER_SEC << " sec." << endl;


	///////////////////  mathcing SNP  ////////////////////////////////////////////////////////////
	CharacterVector rs_12 = intersect(rsname1, rsname2);
	CharacterVector rs_123 = intersect(rs_12, rsname3);
	CharacterVector rs_1235 = intersect(rs_123, rsname5);
	//cout << "rsname1 " << rsname1.length() << endl;
	//cout << "rsname2 " << rsname2.length() << endl;
	//cout << "rsname3 " << rsname3.length() << endl;
	//cout << "rsname5 " << rsname5.length() << endl;
	//cout << "rs_12 " << rs_12.length() << endl;
	//cout << "rs_123 " << rs_123.length() << endl;
	//cout << "rs_1235 " << rs_1235.length() << endl;


	//////////////////  rsname_4use ///////////////////////////////////////////////////////////////
	IntegerVector idxin = match(rsname3, rs_1235);
	CharacterVector rsname_4use = rsname3[Rcpp::is_na(idxin) == false]; // rsname_4use in the order of panel.
	IntegerVector bp3_4use = bp3[Rcpp::is_na(idxin) == false]; // bp3 in the order of panel
	IntegerVector chr3_4use = chr3[Rcpp::is_na(idxin) == false]; // chr3 in the order of panel

	//----------------------------------------------------------------------------------------------//
	printf("Start matching SNPs with summary stat for eqtls. \n");
	IntegerVector idxin1 = match(rsname_4use, rsname1);  //index for SNPs in eqtls

	//----------------------------------------------------------------------------------------------//
	printf("Start matching SNPs with summary stat for GWAS. \n");
	IntegerVector idxin2 = match(rsname_4use, rsname2);  //index for SNPs in GWAS

	//----------------------------------------------------------------------------------------------//
	printf("Start matching SNPs with first reference panel. \n");
	IntegerVector idxin3 = match(rsname_4use, rsname3); //index for SNPs in panel SNPs

	//----------------------------------------------------------------------------------------------//
	printf("Start matching SNPs with second reference panel. \n");
	IntegerVector idxin5 = match(rsname_4use, rsname5); //index for SNPs in panel SNPs

	cout << "Size of matched SNPs: " << rsname_4use.size() << endl; //compare keepIndx in R: Height_1000Q_match.R


	//////////////////   correct direction of minor allele  ////////////////////////////////////////////
	uvec tmp;

	uvec idx1 = as<uvec>(idxin1) -1;
	tmp = as<uvec>(A11);
	uvec A11_ = tmp.elem(idx1);
	tmp = as<uvec>(A12);
	uvec A12_ = tmp.elem(idx1);

	uvec idx2 = as<uvec>(idxin2) -1;
	tmp = as<uvec>(A21);
	uvec A21_ = tmp.elem(idx2);
	tmp = as<uvec>(A22);
	uvec A22_ = tmp.elem(idx2);

	uvec idx3 = as<uvec>(idxin3) -1;
	tmp = as<uvec>(A31);
	uvec A31_ = tmp.elem(idx3);
	tmp = as<uvec>(A32);
	uvec A32_ = tmp.elem(idx3);

	uvec idx5 = as<uvec>(idxin5) -1;
	tmp = as<uvec>(A51);
	uvec A51_ = tmp.elem(idx5);
	tmp = as<uvec>(A52);
	uvec A52_ = tmp.elem(idx5);

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
	idx = find(A51_ > 85);
	A51_.elem(idx) = A51_.elem(idx) - 32;
	idx = find(A52_ > 85);
	A52_.elem(idx) = A52_.elem(idx) - 32;


	//A11_: replace T with A,replace G with C
	T2AG2C(A11_);

	//A12_: replace T with A,replace G with C
	T2AG2C(A12_);

	//A21_: replace T with A,replace G with C
	T2AG2C(A21_);

	//A22_: replace T with A,replace G with C
	T2AG2C(A22_);

	//A31_: replace T with A, replace G with C
	T2AG2C(A31_);

	//A32_: replace T with A,replace G with C
	T2AG2C(A32_);

	//A51_: replace T with A,replace G with C
	T2AG2C(A51_);

	//A52_: replace T with A,replace G with C
	T2AG2C(A52_);

	//remaining index having the same allele
	uvec idxtmp1, idxtmp2, idxtmp5;
	idxtmp1 = find((A31_ + A32_) == (A11_ + A12_));
	idxtmp2 = find((A31_ + A32_) == (A21_ + A22_));
	idxtmp5 = find((A31_ + A32_) == (A51_ + A52_));
	//cout << A11_.n_elem << "  " << A21_.n_elem << "  " << A31_.n_elem << "  " << A51_.n_elem << endl;
	//cout << idxtmp1.n_elem << "   " << idxtmp2.n_elem << "   " << idxtmp5.n_elem << endl;
	uvec idx4_tmp = intersect(idxtmp1, idxtmp2);
	uvec idx4 = intersect(idx4_tmp, idxtmp5);

	uvec A11_r = A11_.elem(idx4), A12_r = A12_.elem(idx4);
	uvec A21_r = A21_.elem(idx4), A22_r = A22_.elem(idx4);
	uvec A31_r = A31_.elem(idx4), A32_r = A32_.elem(idx4);
	uvec A51_r = A51_.elem(idx4), A52_r = A52_.elem(idx4);

	IntegerVector idxtmp = wrap(idx4);
	CharacterVector rsname_4use_r = rsname_4use[idxtmp];
	IntegerVector bp_4use_r = bp3_4use[idxtmp];
	IntegerVector chr_4use_r = chr3_4use[idxtmp];

	cout << "Size of A11_r is " << A11_r.n_elem << endl;
	cout << "Size of A21_r is " << A21_r.n_elem << endl;
	cout << "Size of A31_r is " << A31_r.n_elem << endl;
	cout << "Size of A51_r is " << A51_r.n_elem << endl;

	////////////////////////////////////////////////////////////////////////////////////////////////
	IntegerVector idxinFile2_ = match(rsname_4use_r, rsname2) - 1;  //index for SNPs in GWAS
	IntegerVector idxinFile3_ = match(rsname_4use_r, rsname3) - 1;  //index for SNPs in first reference panel
	IntegerVector idxinFile5_ = match(rsname_4use_r, rsname5) - 1;  //index for SNPs in second reference panel

	uvec idxinFile2 = as<uvec>(idxinFile2_);
	uvec idxinFile3 = as<uvec>(idxinFile3_);
	uvec idxinFile5 = as<uvec>(idxinFile5_);

	uvec idx21 = find(A11_r != A21_r);
	cout << "Size of SNPs with opposite allele direction for GWAS compare with eQTL " << idx21.n_elem << endl;
	uvec tmp1 = idx2(idx4);
	GWAS_sum(tmp1(idx21), zeros<uvec>(1)) = -GWAS_sum(tmp1(idx21), zeros<uvec>(1));

	uvec idx22 = find(A11_r != A31_r);
	cout << "Size of SNPs with opposite allele direction for First reference panel compare with eQTL " << idx22.n_elem << endl;
	uvec tmp2 = idx3(idx4);
	X3_point->cols(tmp2(idx22)) = 2 - X3_point->cols(tmp2(idx22));

	uvec idx25 = find(A11_r != A51_r);
	cout << "Size of SNPs with opposite allele direction for Second reference panel compare with eQTL " << idx25.n_elem << endl;
	uvec tmp5 = idx5(idx4);
	X5_point->cols(tmp5(idx25)) = 2 - X5_point->cols(tmp5(idx25));


	///////////////////  begin to fit the model //////////////////////////
	mat X = StandardX(X3_point);
	mat X2 = StandardX(X5_point);

	// set the parameter
	Options_CoMM_S4* lp_opt = new Options_CoMM_S4(10000, 0, 10, 1e-5, 0);
	Options_CoMM_S4* lp_opt1 = new Options_CoMM_S4(10000, 0, 10, 1e-5, 1);

	// number of gene
	int Ngene = gene_index.n_rows;
	vec stat = zeros(Ngene, 1);
	vec alphag = zeros(Ngene, 1);
	uvec NumSNP = zeros<uvec>(Ngene, 1);
	field<uvec> idx_all(Ngene);
	field<uvec> idx31_all(Ngene);

	t1 = clock();
	cout << "### Precompute the active genes ... ";
	uvec idx_active_gene = zeros<uvec>(0);
	uvec g_tmp(1);
	for (uword g = 0; g < Ngene; g++){

		int Nsnp_gene_eQTL = gene_index(g, 1) - gene_index(g, 0) + 1;
		IntegerVector idx_chr = wrap(linspace<uvec>(gene_index(g, 0), gene_index(g, 1), Nsnp_gene_eQTL) - 1);
		CharacterVector snp_gene_eqtls = eqtls_snp[idx_chr];
		IntegerVector idx_4use_r = match(snp_gene_eqtls, rsname_4use_r);
		CharacterVector snp_gene = snp_gene_eqtls[Rcpp::is_na(idx_4use_r) == false]; //snp in the order of snp_each_gene
		int Nsnp_gene = snp_gene.length();

		if (Nsnp_gene != 0){
			if (Nsnp_gene > 1){
				g_tmp(0) = g;
				idx_active_gene = join_cols(idx_active_gene, g_tmp);
			}
		}

		IntegerVector idx_4use_r_order_ = match(snp_gene, rsname_4use_r);
		uvec idx_4use_r_order = as<uvec>(idx_4use_r_order_) -1;

		idx_all(g) = idx_4use_r_order;

		IntegerVector index1 = match(snp_gene, snp_gene_eqtls);
		uvec idx31 = as<uvec>(index1) -1;

		idx31_all(g) = idx31;

	}
	cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;


	CharacterVector gene_type1 = eqtls_gene[as<Rcpp::IntegerVector>(wrap(idx_active_gene))];

	// combine gene info into mat
	uword Ngene_active = idx_active_gene.n_elem;
	mat gene_info = zeros<mat>(Ngene_active, 3);

	// testing 
	cout << "### Start running CoMM_2S2 on active genes ... " << endl;
	t1 = clock();

	mat out_param0 = -99 * ones<mat>(Ngene_active, 3);

	//set parallel structure object
	parGene_CoMM_2S2 parObj(gene_index, idx31_all, idx_all, idxinFile2, idxinFile3, idxinFile5, GWAS_sum, eqtls_sum, X, X2, px, out_param0, (int)Ngene_active, idx_active_gene, lam);

	const int n_thread = coreNum;
	std::vector<std::thread> threads(n_thread);
	for (int i_thread = 0; i_thread < n_thread; i_thread++){
		threads[i_thread] = std::thread(&parGene_CoMM_2S2::update_by_thread_CoMM_2S2, &parObj, i_thread);
	}

	for (int i = 0; i < n_thread; i++){
		threads[i].join();
	}

	cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;

	mat out_param = -99 * ones<mat>(Ngene, 3);
	out_param.rows(idx_active_gene) = parObj.out_param;


	List output = List::create(
		Rcpp::Named("genetype1") = eqtls_gene,
		Rcpp::Named("gene_type1") = gene_type1,
		Rcpp::Named("out_param") = out_param,
		Rcpp::Named("out_param0") = parObj.out_param);
	return output;

}