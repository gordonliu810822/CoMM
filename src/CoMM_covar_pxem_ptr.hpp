#ifndef CoMM_covar_pxem_ptr_hpp
#define CoMM_covar_pxem_ptr_hpp

#include <RcppArmadillo.h>
//#include <Rcpp.h>
//#include <omp.h>
#include <math.h>


using namespace Rcpp;
using namespace arma;
using namespace std;

void CoMM_covar_pxem_ptr(const arma::vec& y, const arma::vec& z, const arma::mat& x1,  const arma::mat& x2,  const arma::mat& w1,  const arma::mat& w2, double& sigma2beta,
                         double& sigma2y, double& sigma2z, arma::vec& alpha0, double& alpha, double& gam, arma::vec& beta0, double& loglik_max, int& iteration,
                         const int& constr, const double& epsStopLogLik, const int& maxIter);


#endif /* CoMM_covar_pxem_ptr_hpp */
