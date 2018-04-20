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
#include "lmm_covar_pxem.hpp"
#include "AUDI_covar_pxem.hpp"
#include "readExprFile.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;

Rcpp::List getColNum_Header(std::string filename, char delimiter);

Rcpp::List dataLoader(std::string stringname1, std::string stringname2, std::string stringname3, 
	std::string stringname4, std::string stringname5, int whCol);

#endif /* data_loader_hpp */
