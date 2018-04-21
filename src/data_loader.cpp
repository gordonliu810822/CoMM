#include <RcppArmadillo.h>
//#include <Rcpp.h>
//#include <omp.h>
#include <stdio.h>
#include <bitset>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <armadillo_bits/config.hpp>
#include "plinkfun.hpp"
#include "lmm_covar_pxem.hpp"
#include "CoMM_covar_pxem.hpp"
#include "readExprFile.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH, bigmemory)]]
// [[Rcpp::plugins(cpp11)]]

Rcpp::List getColNum_Header(std::string filename, char delimiter){
	//cout << "Count columns in the file: ";

	std::ifstream myfile (filename.c_str());
	std::string line, temp;
	
	if (myfile.is_open()){
		getline(myfile, line);
	}

	stringstream ss(line); // Turn the string into a stream.
	string tok;
	CharacterVector header;

	while(getline(ss, tok, delimiter)) {
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
	}
	while(iss);

	//cout << columns << endl;
	
	List output = List::create(//Rcpp::Named("NumGene") = NumGene,
		Rcpp::Named("columns") = columns,
		Rcpp::Named("header") = header);

	return output;
}

CharacterVector charv_subset2(CharacterVector x, uvec idx){
	CharacterVector v(idx.n_elem);
	for (unsigned int i = 0; i < idx.n_elem; i++){
		v(i) = x(idx(i));
	}
	return v;
}

// format of covariate file is "FID IID, covar1, covar2, covar3 ..."
// [[Rcpp::export]]
Rcpp::List getCovarFile(std::string filename, char delimiter, int ncols, int nrows){
  std::ifstream myfile (filename.c_str());
  mat covar(nrows, ncols-2);
  std::string line;
  CharacterVector IID(nrows), FID(nrows);

  //clock_t t1 = clock();

  int nrow_ind = 0;
  vector <string> tmp;

  if (myfile.is_open()){
	  while (nrow_ind < nrows ){

		  getline(myfile, line);

		  boost::split(tmp, line, boost::is_any_of(" \t *"));

		  IID(nrow_ind) = tmp[1];
		  FID(nrow_ind) = tmp[0];

		  //   cout << ncols-ncols_omit << ";" << nrow_ind << endl;


		  for (int j = 0; j < ncols - 2; j++){
			  if (tmp[j + 2].compare("NA") == 0){
				  covar(nrow_ind, j) = log(-10);
			  }
			  else {
				  covar(nrow_ind, j) = atof(tmp[j + 2].c_str());
			  }
		  }

		  nrow_ind++;

	  }
  }

  List output = List::create(Rcpp::Named("IID") = IID,
	  Rcpp::Named("FID") = FID,
      Rcpp::Named("covar") = covar);
  return output;
}



// [[Rcpp::export]]
Rcpp::List dataLoader(std::string stringname1, std::string stringname2, std::string stringname3, 
	std::string stringname4, std::string stringname5, int whCol){
	//match SNPs in file 1 and file 2 GWAS (common SNPs in x1 and x2 in columns)
	// plink file 1: stringname1; plink file 2: stringname2; expression file: stringname3
	// covariates file for file 1: stringname4; covariates file for file 2: stringname5
	// whCol: which pheno is used (1 stands for first pheno (6th column of fam file and so on )
	cout << "## Start matching SNPs in plink files (1, 2). " << endl;

	List tmp =  match_SNPs(stringname1, stringname2);
	uvec idxinFile1 = tmp["idxinFile1"];
	uvec idxinFile2 = tmp["idxinFile2"];
	vec ind = tmp["indicator"];
	CharacterVector rsname_4use_r = tmp["rsname_4use_r"];
	uvec chr_4use_r = tmp["chr_4use_r"];
	uvec bp_4use_r = tmp["bp_4use_r"];
	uvec A1_1_r = tmp["A1_1_r"], A1_2_r = tmp["A1_2_r"], A2_1_r = tmp["A2_1_r"], A2_2_r = tmp["A2_2_r"];

	// load IID in file 1 (fam file 1)
	cout << "## Start loading fam files (1, 2). " << endl;
	string famfile1 = stringname1;
	famfile1 += ".fam";
	int N1 = getLineNum(famfile1);

	IntegerVector tmp1(N1);
	CharacterVector FID_1(N1), IID_1(N1);
	NumericVector tmp2(N1);
	ReadPlinkFamFile(famfile1, FID_1, IID_1, tmp1, tmp2, N1);

	// load pheno in file 2 (fam file 2)
	string famfile2 = stringname2;
	famfile2 += ".fam";
	int N2 = getLineNum(famfile2);

	IntegerVector sex_2(N2);
	NumericVector pheno_2(N2);
	CharacterVector FID_2(N2), IID_2(N2);
	// ReadPlinkFamFile(famfile2, FID_2, IID_2, sex_2, pheno_2, N2);
	ReadPlinkFamFile2(famfile2, FID_2, IID_2, pheno_2, N2, whCol);

	vec y = as<vec>(pheno_2);

	// read expression file
	cout << "## Start loading expression file. " << endl;
	char delimiter = '\t';
	tmp = getColNum_Header(stringname3, delimiter);
	int Nindiv = tmp["columns"];
	CharacterVector header = tmp["header"];
	int Ngene = getLineNum(stringname3) - 1; // header is first line;
	
	// IID from header of expr file
	uvec idxtmp = linspace<uvec>(6, Nindiv-1, Nindiv - 6);
	CharacterVector IID_expr = charv_subset2(header, idxtmp);
	
	// match file 1 with expression file
	CharacterVector indiv_inter = intersect(IID_1, IID_expr);
	IntegerVector idxint = match(IID_1, indiv_inter); // put the IID in the order of file 1
	CharacterVector indiv_4use = IID_1[Rcpp::is_na(idxint) == false];
	IntegerVector idxin1 = match(indiv_4use, IID_1) -1;  //index for individual in file 1
	IntegerVector idxin2 = match(indiv_4use, IID_expr) -1;  //index for individual in expr file

	//read expr file
	List expr_tmp = getExprFile2(stringname3, delimiter, Nindiv, Ngene + 1, 6);
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
	exprTmp = & expr;
	expr_used = exprTmp -> cols(as<uvec>(idxin2));

	clock_t t1 = clock();
	cout << "## Start loading genotype files 1, " ;
	//read file 1 (plink bim bed fam files)
	string bimfile1 = stringname1;
	bimfile1 += ".bim";
	long int P1 =  getLineNum(bimfile1);
	cout << N1 << "-by-" << P1;
	//unsigned* X1tmp = new unsigned[ N1 * P1]; 
	arma::Mat<unsigned> X1(N1,P1);
	//cout << " break 1: " << N1 << "X" << P1 << endl;
	readPlink(stringname1,N1, P1, X1.memptr());
	
	// subset overlapped SNPs and individuals
	//arma::Mat<unsigned>* Xdata = new arma::Mat<unsigned>(X1tmp, N1, P1, false, false);
	/*arma::Mat<unsigned>* Xdata = new arma::Mat<unsigned>(X1tmp, N1, P1, false, false);
	//arma::Mat<unsigned>* Xmtmp;
	arma::Mat<unsigned> subX, X1;

	subX = Xdata->cols(idxinFile1 - 1);
	Xdata = &subX;
	Xdata->replace(3, 0);
	X1 = Xdata->rows(as<uvec>(idxin1));*/

	// arma::Mat<unsigned> X1(X1tmp, N1, P1, false, false);
	//cout << " break 2: " << X1.n_rows << "X" << X1.n_cols << endl;

	X1 = X1.cols(idxinFile1 - 1);
	X1 = X1.rows(as<uvec>(idxin1));
	// replace NA with 0 dosage
	X1.replace(3, 0);

	//delete[] X1tmp;
	cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
	
	t1 = clock();
	cout << "## Start loading genotype files 2, ";
	//read file 2 (plink bim bed fam files)
	//cout << endl;
	string bimfile2 = stringname2;
	bimfile2 += ".bim";
	//cout << "break 1... " << endl;
	long int P2 =  getLineNum(bimfile2);
	cout << N2 << "-by-" << P2;
	//unsigned* X2tmp = new unsigned[ N2 * P2]; 
	arma::Mat<unsigned> X2(N2,P2);
	//cout << "break 3..." << endl;
	readPlink(stringname2,N2, P2, X2.memptr());

	// subset overlapped SNPs and individuals
    /*arma::Mat<unsigned>* Xdata2 = new arma::Mat<unsigned>(X2tmp, N2, P2, false,false);
	Xdata2->replace(3, 0);
	X2 = Xdata2 ->cols(idxinFile2-1);*/
	// arma::Mat<unsigned> X2(X2tmp, N2, P2, false, false);
	X2 = X2.cols(idxinFile2 -1);
	// replace NA with 0 dosage
	X2.replace(3, 0);

	//delete[] X2tmp;

	//change the reference allele in file 2
	idxtmp = linspace<uvec>(1, idxinFile1.size(), idxinFile1.size());
	uvec idx = find(ind == -1);
	uvec idxc = idxtmp.elem(idx);

	X2.cols(idxc-1) = 2 - X2.cols(idxc-1);
	cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
	//int Nindiv = getColNum(stringname3);

	cout << "## Start loading covariates files ... " << endl;
	// load covariates file w2
	mat covar2;
	if (!stringname5.empty()){
		tmp = getColNum_Header(stringname5, delimiter);
		int Ncovar = tmp["columns"];
		tmp = getCovarFile(stringname5, delimiter, Ncovar, N2);
		mat covartmp = tmp["covar"];
		
		mat w2one = ones<mat>(N2, 1);
		//cout << endl << "break ..." << covartmp.n_rows << ";" << covartmp.n_cols << ";" << N2 << endl;

		covar2 = join_rows(w2one, covartmp);
	}
	else {
		covar2 = ones<mat>(N2,1);
	}
	
	// load covariates file w1
	mat covar1;
	//cout << "break ..." << covar2.n_rows << ";" << covar2.n_cols << ";" << N2 << ";" << !stringname5.empty() << endl;
	CharacterVector IID_w1;
	if (!stringname4.empty()){
		tmp = getColNum_Header(stringname4, delimiter);
		int Ncovar = tmp["columns"];
		tmp = getCovarFile(stringname4, delimiter, Ncovar, N1);
		mat covar1 = tmp["covar"];
		
		IID_w1 = tmp["IID"];

		covar1 = join_rows(ones<mat>(X1.n_rows, 1), covar1.rows(as<uvec>(idxin1)));
	}
	else {
		covar1 = ones<mat>(X1.n_rows,1);
	}
	cout << "## End loading files ... " << endl;
	List output = List::create(/*Rcpp::Named("Numgene") = Ngene,
		Rcpp::Named("Nindiv") = Nindiv,
		Rcpp::Named("IID_1") = IID_1,
		Rcpp::Named("IID_expr") = IID_expr,
		Rcpp::Named("idxinFile1") = idxinFile1,
		Rcpp::Named("idxinFile2") = idxinFile2,
		Rcpp::Named("idxin1") = idxin1,
		Rcpp::Named("idxin2") = idxin2,
		Rcpp::Named("A1_1_r") = A1_1_r,
		Rcpp::Named("A1_2_r") = A1_2_r,
		Rcpp::Named("A2_1_r") = A2_1_r,
		Rcpp::Named("A2_2_r") = A2_2_r,
		Rcpp::Named("indicator") = ind,
		Rcpp::Named("expr") = expr,*/
		//Rcpp::Named("idxc") = idxc,
		//Rcpp::Named("indiv_4use") = indiv_4use,
		//Rcpp::Named("IID_w1") = IID_w1,
		//Rcpp::Named("idxin1") = idxin1,
		Rcpp::Named("X1") = X1,
		Rcpp::Named("X2") = X2,
		Rcpp::Named("y") = y,
		Rcpp::Named("covar1") = covar1,
		Rcpp::Named("covar2") = covar2,
		Rcpp::Named("expr_used") = expr_used,
		Rcpp::Named("rsname_4use_r") = rsname_4use_r,
		Rcpp::Named("chr_4use_r") = chr_4use_r,
		Rcpp::Named("bp_4use_r") = bp_4use_r,
		Rcpp::Named("targetID") = targetID,
		Rcpp::Named("genetype1") = genetype1,
		Rcpp::Named("genetype2") = genetype2,
		Rcpp::Named("lower") = lower,
		Rcpp::Named("upper") = upper,
		Rcpp::Named("chr_expr") = chr_expr);
		//Rcpp::Named("test") = ((&X1) ->cols(conv_to<uvec>::from(0))));

	return output;
}
