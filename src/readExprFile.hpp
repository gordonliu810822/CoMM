#ifndef readExprFile_hpp
#define readExprFile_hpp

#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include "plinkfun.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace Rcpp;
using namespace arma;

Rcpp::List getExprFile2(std::string filename, char delimiter, int ncols, int nrows, int ncols_omit);

#endif /* readExprFile_hpp */
