#ifndef lmm_covar_pxem_hpp
#define lmm_covar_pxem_hpp

#include <RcppArmadillo.h>
//#include <armadillo>
//#include <Rcpp.h>
//#include <omp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

List lmm_pxem(arma::vec y, arma::mat w, arma::mat x, int maxIter);

#endif /* lmm_covar_pxem_hpp */
