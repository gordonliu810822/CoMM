#ifndef data_loader_hpp
#define data_loader_hpp


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

#include "plinkfun.hpp"
//#include "lmm_covar_pxem.hpp"
#include "lmm_covar_pxem_ptr.hpp"
//#include "AUDI_covar_pxem.hpp"
#include "CoMM_covar_pxem_ptr.hpp"
#include "readExprFile.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;

Rcpp::List getColNum_Header(std::string filename, char delimiter);
Rcpp::List getCovarFile(std::string filename, char delimiter, int ncols, int nrows);
Rcpp::List dataLoader1(std::string stringname1, std::string stringname2, std::string stringname3, int whCol);
Rcpp::List dataLoader2(std::string stringname1, std::string stringname2, std::string stringname4, std::string stringname5,
                       Mat<unsigned>& X1, Mat<unsigned>& X2, uvec& idxinFile1, uvec& idxin1, uvec& idxinFile2, uvec& idxin2, vec& ind);
//Rcpp::List dataLoader3(std::string stringname1, std::string stringname2, std::string stringname4, std::string stringname5, char* X1, char* X2)；

#endif /* data_loader_hpp */
