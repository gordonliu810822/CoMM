#ifndef CoMM_S2_testing_hpp
#define CoMM_S2_testing_hpp

#include <RcppArmadillo.h>
//#include <Rcpp.h>
//#include <omp.h>
#include <math.h>


using namespace Rcpp;
using namespace arma;
using namespace std;

List read_GWAS(string filename, int P);

Rcpp::List ReadSNPinfo(std::string stringname, IntegerVector A1, IntegerVector A2, CharacterVector rsname,
	IntegerVector chr, IntegerVector bp, NumericVector morgan, int N);

void T2AG2C(uvec& A31_);

Rcpp::List getColNum_Header2(std::string filename, char delimiter);

CharacterVector charv_subset2_(CharacterVector x, uvec idx);

void ReadPlinkFamFile2(std::string stringname, CharacterVector FID, CharacterVector IID, IntegerVector sex,
	NumericVector pheno, int N);

#endif /* CoMM_covar_pxem_ptr_hpp */
