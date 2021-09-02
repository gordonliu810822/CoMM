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
#include "lmm_covar_pxem_ptr.hpp"
#include "CoMM_covar_pxem_ptr.hpp"
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

List match_SNPs(std::string stringname1, std::string stringname2){
    
    // stringname1: prefix for plink file 1; start working on plink file 1 (1000 G with expression data)
    string bimfile1 = stringname1;
    bimfile1 += ".bim";
    //cout << "break 1: " << bimfile1  << endl;
    int P1 = getLineNum(bimfile1);
    cout << "## Number of SNPs (plink file 1):" << P1 << endl;
    
    IntegerVector A1_1(P1), A2_1(P1);
    CharacterVector rsname_1(P1);
    IntegerVector chr_1(P1), bp_1(P1);
    NumericVector morgan_1(P1);
    
    ReadPlinkBimFile(bimfile1, A1_1, A2_1, rsname_1, chr_1, bp_1, morgan_1, P1);
    // cout << rsname_1(0) << ";" << A1_1(0) << ";" << A2_1(0) << endl;
    
    // stringname2: prefix for plink file 2; start working on plink file 2 (GWAS data with trait)
    string bimfile2 = stringname2;
    bimfile2 += ".bim";
    
    int P2 = getLineNum(bimfile2);
    cout << "## Number of SNPs (plink file 2):" << P2 << endl;
    
    IntegerVector A1_2(P2), A2_2(P2);
    CharacterVector rsname_2(P2);
    IntegerVector chr_2(P2), bp_2(P2);
    NumericVector morgan_2(P2);
    
    ReadPlinkBimFile(bimfile2, A1_2, A2_2, rsname_2, chr_2, bp_2, morgan_2, P2);
    // cout << rsname_2(0) << ";" << chr_2(0) << ";" << bp_2(0) << ";" << A1_2(0) << ";" << A2_2(0) << endl;
    
    // mathcing panel SNPs in file 1 and file 2 with correction for direction of ref allele
    // rsname in both files in the order of the first file.
    CharacterVector rs_inter = intersect(rsname_1, rsname_2);
    IntegerVector idxin = match(rsname_1, rs_inter); //index for SNPs in file 1
    CharacterVector rsname_4use = rsname_1[Rcpp::is_na(idxin) == false];
    IntegerVector chr_4use = chr_1[Rcpp::is_na(idxin) == false];
    IntegerVector bp_4use = bp_1[Rcpp::is_na(idxin) == false];
    
    // match snps (rsname_4use; rsname_2: file 2)
    IntegerVector idxin2 = match(rsname_4use, rsname_2);  //index for SNPs in file 2
    // match snps (rsname_4use; rsname_1: file 1)
    IntegerVector idxin1 = match(rsname_4use, rsname_1);  //index for SNPs in file 1
    
    /*IntegerVector chr_4use_tmp = chr_2[idxin2];
     IntegerVector bp_4use_tmp = bp_2[idxin2];
     
     vec idxtmp = as<vec>(chr_4use) - as<vec>(chr_4use_tmp);
     cout << "check the quality: " << sum(idxtmp) << endl;
     
     cout << "Size of matched SNPs: " << rsname_4use.size() << endl;*/
    
    // convert ascii letter to uvec and work on overlapped SNPs
    uvec idx = as<uvec>(idxin1) -1;
    uvec tmp1 = as<uvec>(A1_1);
    uvec A1_1_ = tmp1.elem(idx);
    tmp1 = as<uvec>(A2_1);
    uvec A2_1_ = tmp1.elem(idx);
    
    idx = as<uvec>(idxin2) -1;
    uvec tmp2 = as<uvec>(A1_2);
    uvec A1_2_ = tmp2.elem(idx);
    tmp2 = as<uvec>(A2_2);
    uvec A2_2_ = tmp2.elem(idx);
    
    //ascii: (A:65; C:67; G:71; T:84) (a:97;c:99,g:103;t:116)
    /*//replace lower letter to upper letter in A1_1_ and A2_1_;
     idx = find(A1_1_ > 85);
     A1_1_.elem(idx) = A1_1_.elem(idx) - 32;
     idx = find(A2_1_ > 85);
     A2_1_.elem(idx) = A2_1_.elem(idx) - 32;
     
     idx = find(A1_2_ > 85);
     A1_2_.elem(idx) = A1_2_.elem(idx) - 32;
     idx = find(A2_2_ > 85);
     A2_2_.elem(idx) = A2_2_.elem(idx) - 32;*/
    
    //compare A1_1_ A1_2_, A2_1_ A2_2_
    //A1_1_: replace T with A,replace G with C
    idx = find(A1_1_ == 84);
    uvec idxrepl1(idx.n_elem);
    idxrepl1.fill(65);
    A1_1_.elem(idx) = idxrepl1;
    
    idx = find(A1_1_ == 71);
    uvec idxrepl2(idx.n_elem);
    idxrepl2.fill(67);
    A1_1_.elem(idx) = idxrepl2;
    
    //A1_2_: replace T with A,replace G with C
    idx = find(A1_2_ == 84);
    uvec idxrepl3(idx.n_elem);
    idxrepl3.fill(65);
    A1_2_.elem(idx) = idxrepl3;
    
    idx = find(A1_2_ == 71);
    uvec idxrepl4(idx.n_elem);
    idxrepl4.fill(67);
    A1_2_.elem(idx) = idxrepl4;
    
    //A2_1_: replace T with A,replace G with C
    idx = find(A2_1_ == 84);
    uvec idxrepl5(idx.n_elem);
    idxrepl5.fill(65);
    A2_1_.elem(idx) = idxrepl5;
    
    idx = find(A2_1_ == 71);
    uvec idxrepl6(idx.n_elem);
    idxrepl6.fill(67);
    A2_1_.elem(idx) = idxrepl6;
    
    //A2_2_: replace T with A,replace G with C
    idx = find(A2_2_ == 84);
    uvec idxrepl7(idx.n_elem);
    idxrepl7.fill(65);
    A2_2_.elem(idx) = idxrepl7;
    
    idx = find(A2_2_ == 71);
    uvec idxrepl8(idx.n_elem);
    idxrepl8.fill(67);
    A2_2_.elem(idx) = idxrepl8;
    
    // remove index;
    idx = find((A1_1_ + A2_1_) == (A1_2_ + A2_2_));
    uvec tmp3 = as<uvec>(idxin2);
    uvec idxin22 = tmp3.elem(idx);
    tmp3 = as<uvec>(idxin1);
    uvec idxin11 = tmp3.elem(idx);
    
    uvec A1_1_r = A1_1_.elem(idx), A2_1_r = A2_1_.elem(idx);
    uvec A1_2_r = A1_2_.elem(idx), A2_2_r = A2_2_.elem(idx);
    
    uvec chr_tmp = as<uvec>(chr_4use);
    uvec chr_4use_r = chr_tmp.elem(idx);
    uvec bp_tmp = as<uvec>(bp_4use);
    uvec bp_4use_r = bp_tmp.elem(idx);
    
    CharacterVector rsname_4use_r = charv_subset(rsname_4use, idx);
    
    vec ind(A1_1_r.n_elem);
    idx = find(A1_1_r == A1_2_r);
    ind.ones();
    ind = -ind;
    ind.elem(idx).ones();
    
    cout << "Number of matched SNPs (plink file 1 and 2):" << ind.size() << endl;
    //cout << "direction (compare with trait_1000Q_match.R): " << sum(ind) << endl; //compare sum(ind) in R: Height_1000Q_match.R
    //cout << "Size of matched SNPs (remove ambiguous SNPs): " << A1_1_r.n_elem << endl;
    
    List output = List::create(Rcpp::Named("chr_4use_r") = chr_4use_r,
                               Rcpp::Named("A1_1_r") = A1_1_r,
                               Rcpp::Named("A1_2_r") = A1_2_r,
                               Rcpp::Named("A2_1_r") = A2_1_r,
                               Rcpp::Named("A2_2_r") = A2_2_r,
                               Rcpp::Named("bp_4use_r") = bp_4use_r,
                               Rcpp::Named("rsname_4use_r") = rsname_4use_r,
                               Rcpp::Named("indicator") = ind,
                               Rcpp::Named("idxinFile1") = idxin11,
                               Rcpp::Named("idxinFile2") = idxin22);
    
    return output;
}

Rcpp::List dataLoader1(std::string stringname1, std::string stringname2, std::string stringname3, int whCol){
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
    
    List output = List::create(Rcpp::Named("y") = y,
                               Rcpp::Named("expr_used") = expr_used,
                               Rcpp::Named("rsname_4use_r") = rsname_4use_r,
                               Rcpp::Named("chr_4use_r") = chr_4use_r,
                               Rcpp::Named("bp_4use_r") = bp_4use_r,
                               Rcpp::Named("targetID") = targetID,
                               Rcpp::Named("genetype1") = genetype1,
                               Rcpp::Named("genetype2") = genetype2,
                               Rcpp::Named("lower") = lower,
                               Rcpp::Named("upper") = upper,
                               Rcpp::Named("chr_expr") = chr_expr,
                               Rcpp::Named("idxin1") = as<uvec>(idxin1),
                               Rcpp::Named("idxin2") = as<uvec>(idxin2),
                               Rcpp::Named("idxinFile1") = idxinFile1,
                               Rcpp::Named("idxinFile2") = idxinFile2,
                               Rcpp::Named("ind") = ind,
							   Rcpp::Named("indiv_4use") = indiv_4use);
    
    return output;
    
}

List dataLoader2(std::string stringname1, std::string stringname2, std::string stringname4, std::string stringname5, 
				 Mat<unsigned>& X1, Mat<unsigned>& X2, uvec& idxinFile1, uvec& idxin1, uvec& idxinFile2, uvec& idxin2, vec& ind)
{
    //match SNPs in file 1 and file 2 GWAS (common SNPs in x1 and x2 in columns)
    // plink file 1: stringname1; plink file 2: stringname2;
    // covariates file for file 1: stringname4; covariates file for file 2: stringname5
    char delimiter = '\t';
    // load the size of plink file
    string famfile1 = stringname1;
    famfile1 += ".fam";
    int N1 = getLineNum(famfile1);
    string bimfile1 = stringname1;
    bimfile1 += ".bim";
    long int P1 =  getLineNum(bimfile1);
    string famfile2 = stringname2;
    famfile2 += ".fam";
    int N2 = getLineNum(famfile2);
    string bimfile2 = stringname2;
    bimfile2 += ".bim";
    long int P2 =  getLineNum(bimfile2);

    clock_t t1 = clock();
    cout << "## Start loading genotype files 1, " ;
    readPlink(stringname1,N1, P1, X1.memptr());
    
    X1 = X1.cols(idxinFile1 - 1);
    X1 = X1.rows(idxin1);
    // replace NA with 0 dosage
    X1.replace(3, 0);
    
    //delete[] X1tmp;
    cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
    
    t1 = clock();
    cout << "## Start loading genotype files 2, ";
    readPlink(stringname2,N2, P2, X2.memptr());
    
    X2 = X2.cols(idxinFile2 -1);
    // replace NA with 0 dosage
    X2.replace(3, 0);
    
    //delete[] X2tmp;
    
    //change the reference allele in file 2
    uvec idxtmp = linspace<uvec>(1, idxinFile1.size(), idxinFile1.size());
    uvec idx = find(ind == -1);
    uvec idxc = idxtmp.elem(idx);
    
    X2.cols(idxc-1) = 2 - X2.cols(idxc-1);
    cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
    //int Nindiv = getColNum(stringname3);
    
    cout << "## Start loading covariates files ... " << endl;
    // load covariates file w2
    mat covar2;
    if (!stringname5.empty()){
        List tmp = getColNum_Header(stringname5, delimiter);
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
        List tmp = getColNum_Header(stringname4, delimiter);
        int Ncovar = tmp["columns"];
        tmp = getCovarFile(stringname4, delimiter, Ncovar, N1);
        mat covar1 = tmp["covar"];
        
        IID_w1 = tmp["IID"];
        
        covar1 = join_rows(ones<mat>(X1.n_rows, 1), covar1.rows(idxin1));
    }
    else {
        covar1 = ones<mat>(X1.n_rows,1);
    }
    cout << "## End loading files ... " << endl;
    
    List output = List::create(
                               Rcpp::Named("covar1") = covar1,
                               Rcpp::Named("covar2") = covar2);
    
    return output;
}

/*
List dataLoader3(std::string stringname1, std::string stringname2, std::string stringname4, std::string stringname5, char* X1, char* X2)
{
    //match SNPs in file 1 and file 2 GWAS (common SNPs in x1 and x2 in columns)
    // plink file 1: stringname1; plink file 2: stringname2;
    // covariates file for file 1: stringname4; covariates file for file 2: stringname5
    char delimiter = '\t';
    // load the size of plink file
    string famfile1 = stringname1;
    famfile1 += ".fam";
    int N1 = getLineNum(famfile1);
    string bimfile1 = stringname1;
    bimfile1 += ".bim";
    long int P1 =  getLineNum(bimfile1);
    string famfile2 = stringname2;
    famfile2 += ".fam";
    int N2 = getLineNum(famfile2);
    string bimfile2 = stringname2;
    bimfile2 += ".bim";
    long int P2 =  getLineNum(bimfile2);
    
    clock_t t1 = clock();
    cout << "## Start loading genotype files 1, " ;
    readPlink(stringname1,N1, P1, X1);
    
    //X1 = X1.cols(idxinFile1 - 1);
    //X1 = X1.rows(idxin1);
    // replace NA with 0 dosage
    //X1.replace(3, 0);
    
    //delete[] X1tmp;
    cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
    
    t1 = clock();
    cout << "## Start loading genotype files 2, ";
    readPlink(stringname2,N2, P2, X2);
    
    //X2 = X2.cols(idxinFile2 -1);
    // replace NA with 0 dosage
    //X2.replace(3, 0);
    
    //delete[] X2tmp;
    
    //change the reference allele in file 2
    //uvec idxtmp = linspace<uvec>(1, idxinFile1.size(), idxinFile1.size());
    //uvec idx = find(ind == -1);
    //uvec idxc = idxtmp.elem(idx);
    
    //X2.cols(idxc-1) = 2 - X2.cols(idxc-1);
    cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
    //int Nindiv = getColNum(stringname3);
    
    cout << "## Start loading covariates files ... " << endl;
    // load covariates file w2
    mat covar2;
    if (!stringname5.empty()){
        List tmp = getColNum_Header(stringname5, delimiter);
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
        List tmp = getColNum_Header(stringname4, delimiter);
        int Ncovar = tmp["columns"];
        tmp = getCovarFile(stringname4, delimiter, Ncovar, N1);
        mat covar1 = tmp["covar"];
        
        IID_w1 = tmp["IID"];
        
        covar1 = join_rows(ones<mat>(X1.n_rows, 1), covar1.rows(idxin1));
    }
    else {
        covar1 = ones<mat>(X1.n_rows,1);
    }
    cout << "## End loading files ... " << endl;
    
    List output = List::create(
                               Rcpp::Named("covar1") = covar1,
                               Rcpp::Named("covar2") = covar2);
    
    return output;
}

*/
