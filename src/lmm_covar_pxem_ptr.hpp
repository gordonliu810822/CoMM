#ifndef lmm_covar_pxem_ptr_hpp
#define lmm_covar_pxem_ptr_hpp

#include <RcppArmadillo.h>
//#include <armadillo>
//#include <Rcpp.h>
//#include <omp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

//List lmm_pxem(const arma::vec& y, const arma::mat& w, const arma::mat& x, const int maxIter);
void lmm_pxem_ptr(const arma::vec& y, const arma::mat& w, const arma::mat& x, const int& maxIter, 
			  double& sigma2y, double& sigma2beta, arma::vec& beta0, double& loglik_max, 
			  int& iteration, arma::mat& Sigb, arma::vec& mub);

#endif /* lmm_covar_pxem_ptr_hpp */
