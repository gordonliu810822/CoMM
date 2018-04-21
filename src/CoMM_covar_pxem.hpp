#ifndef twas_covar_pxem_hpp
#define twas_covar_pxem_hpp

#include <RcppArmadillo.h>
//#include <Rcpp.h>
//#include <omp.h>
#include <math.h>


using namespace Rcpp;
using namespace arma;
using namespace std;

double loglike_twas(mat R, vec res_y, vec res_z, vec mu, double sigma2beta, double sigma2y, double sigma2z, int n1,
	int n2, int p);
double loglike_twas2(mat w1, mat w2, vec alpha0, vec beta0, mat x1x1t, mat x1x2t, mat x2x2t,
	double alpha, double sigma2beta, double sigma2y, double sigma2z, vec y, vec z, int n1,
	int n2, int p);
List CoMM_covar_pxem(arma::vec y, arma::vec z, arma::mat x1, arma::mat x2, arma::mat w1, arma::mat w2, double sigma2beta,
	double sigma2y, arma::vec beta0, int constr, double epsStopLogLik, int maxIter);

#endif /* twas_covar_pxem_hpp */
