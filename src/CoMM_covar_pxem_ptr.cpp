#include <RcppArmadillo.h>
//#include <armadillo>
//#include <Rcpp.h>
//#include <omp.h>
#include <math.h>
#include "lmm_covar_pxem_ptr.hpp"
#include "lmm_pxem.hpp"
using namespace Rcpp;
using namespace arma;
using namespace std;

//const double PI  = 3.141592653589793238463;

void loglike_twas(const arma::mat& R, const arma::vec& res_y, const arma::vec& res_z, const arma::vec& mu, const double& sigma2beta, const double& sigma2y,
					const double& sigma2z, const int& n1, const int& n2, const int& p, double& loglik){
	double Eab;//, loglik;
	Eab = sum(res_y % res_y)/2/sigma2y + sum(res_z % res_z)/2/sigma2z + sum(mu % mu)/2/sigma2beta;

	loglik = -log(sigma2beta)*p*0.5 - log(sigma2y)*n1*0.5 - log(sigma2z)*n2 *0.5 - sum(log(R.diag())) - Eab;

}


void CoMM_covar_pxem_ptr(const arma::vec& y, const arma::vec& z, const arma::mat& x1,  const arma::mat& x2,  const arma::mat& w1,  const arma::mat& w2, double& sigma2beta,
                         double& sigma2y, double& sigma2z, arma::vec& alpha0, double& alpha, double& gam, arma::vec& beta0, double& loglik_max, int& iteration,
                         const int& constr, const double& epsStopLogLik, const int& maxIter){
    int n1 = y.n_elem, n2 = z.n_elem, p1 = x1.n_cols, p2 = x2.n_cols;
    if (p1 != p2){
        perror("The dimensions of x1 and x2 are not matched");
    }
    int p = p1;

    mat x1tx1 = x1.t()*x1, x2tx2 = x2.t()*x2, w1tw1 = w1.t()*w1, w2tw2 = w2.t()*w2, x1tw1 = x1.t()*w1, x2tw2 = x2.t()*w2;
    vec x1ty = x1.t()*y, x2tz = x2.t()*z, w2tz = w2.t()*z, w1ty = w1.t()*y;

    //declaration of variables used within loop
    mat Sigma_inv(p,p), R(p,p), invR(p,p), Sigma(p,p);
    vec mutmp(p), mu(p), res_y(n1), res_z(n2), x1mu(n1), x2mu(n2);
    double ztmp;
    mat I = eye<mat>(p,p);

    //initialization of parameters
    if (beta0.n_elem != w1.n_cols){
        perror("The dimensions in covariates are not matched in w1 and beta0");
    }

    if (alpha0.n_elem != w2.n_cols){
        perror("The dimensions in covariates are not matched in w2 and alpha0");
    }

    alpha0 = solve(w2tw2, w2tz);
    alpha = 0;
    sigma2z = var(z);
    gam = 1;
    double gam2 = gam*gam, alpha2 = alpha*alpha;

    vec loglik(maxIter);// loglik2(maxIter);

    // initialization of likelihood
    loglik(0) = NAN;

    int Iteration = 1;
    for (int iter = 2; iter <= maxIter; iter ++ ) {
        // E-step
        Sigma_inv = x1tx1/sigma2y + x2tx2*(alpha2/sigma2z) + I /sigma2beta;
        //cout << "break : " << Sigma_inv << endl;
        R = chol(Sigma_inv);
        invR = inv(R);
        Sigma = invR*invR.t();
        mutmp = (x1ty - x1tw1*beta0)/sigma2y + (x2tz - x2tw2*alpha0)*(alpha/sigma2z);
        mu = Sigma*mutmp;

        //evaluate incomplete log-likelihood
        x1mu = x1*mu;
        x2mu = x2*mu;
        res_y = y - x1mu*gam - w1*beta0;
        res_z = z - alpha*x2mu - w2*alpha0;


        loglike_twas(R, res_y, res_z, mu, sigma2beta, sigma2y, sigma2z,n1,n2,p, loglik(iter - 1));
        /*loglik2(iter - 1)= loglike_twas2(w1, w2, alpha0, beta0, x1x1t, x1x2t, x2x2t,
         alpha, sigma2beta, sigma2y, sigma2z, y, z, n1, n2, p);*/

        //cout << iter << "-th iter: loglik = " << loglik(iter - 1) << "; diff : " << loglik(iter - 1) - loglik(iter - 2) << endl;
        //printf ("%d-th iter: loglik = %11.8f; diff: %11.10f \n", iter,  loglik(iter - 1), loglik(iter - 1) - loglik(iter - 2));

        if ( loglik(iter - 1) - loglik(iter - 2) < -1e-7){
            perror("The likelihood failed to increase!");
        }
        // M-step
        mat trtmp1 = x1tx1 * Sigma;
        double tr1 = trace(trtmp1);
        gam = sum( (x1ty - x1tw1*beta0) % mu )/(sum(x1mu % x1mu) + tr1);
        gam2 = gam*gam;

        beta0 = solve(w1tw1, w1ty - gam*x1tw1.t()*mu);
        alpha0= solve(w2tw2, w2tz - alpha*(x2tw2.t()*mu));

        ztmp = sum( (z - w2*alpha0)% x2mu);
        mat trtmp2 = x2tx2 * Sigma;
        double tr2 = trace(trtmp2);
        if (constr == 1){
            alpha = 0;
        }
        else {
            alpha = ztmp/(sum(x2mu % x2mu) + tr2);
        }
        alpha2 = alpha*alpha;

        res_y = y - x1mu*gam - w1*beta0;
        res_z = z - alpha*x2mu - w2*alpha0;

        sigma2y = (sum(res_y % res_y) + gam2 * tr1)/n1;
        sigma2z = (sum(res_z % res_z) + alpha2 * tr2)/n2;
        sigma2beta = (sum(mu % mu) + trace(Sigma))/p;

        // Reduction-step
        sigma2beta = gam2 * sigma2beta;
        alpha = alpha / gam;
        alpha2 = alpha*alpha;
        gam = 1;
        gam2 = 1;

        Iteration = iter;
        if (iter > 2){
            if (abs(loglik(iter - 1) - loglik(iter - 2)) < epsStopLogLik) {

                break;
            }
        }
    }

    vec loglik_out;
    int to = Iteration -1;
    loglik_out = loglik.subvec(0, to);

    loglik_max = loglik(to);
    iteration = Iteration -1;

}

//' @author Jin Liu, \email{jin.liu@duke-nus.edu.sg}
//' @title
//' CoMM
//' @description
//' CoMM to dissecting genetic contributions to complex traits by leveraging regulatory information.
//'
//' @param y  gene expression vector.
//' @param z  trait vector.
//' @param x1  normalized genotype (cis-SNPs) matrix for eQTL.
//' @param x2  normalized genotype (cis-SNPs) matrix for GWAS.
//' @param w1  covariates file for eQTL data.
//' @param w2  covariates file for GWAS data, e.g. top 10 PCs.
//' @param constr  indicator for constraint (alpha = 0 if constr = 1).
//' @param epsStopLogLik  convergence criteria (default is 1e-5).
//' @param maxIter  maximum iteration (default is 1000).
//' @param pxem_indicator  indicator for using PX-EM (default is 1). pxem_indicator = 1 for using PX-EM and EM otherwise.
//'
//' @return List of model parameters
//'
//' @examples
//' L = 1;
//' M = 100;
//' rho =0.5;
//' n1 = 350;
//' n2 = 5000;
//'
//' X <- matrix(rnorm((n1+n2)*M),nrow=n1+n2,ncol=M);
//' beta_prop = 0.2;
//' b = numeric(); m = M * beta_prop;
//' b[sample(M,m)] = rnorm(m); h2y = 0.05;
//' y0 <- X%*%b + 6;
//' y  <- y0 + (as.vector(var(y0)*(1-h2y)/h2y))^0.5*rnorm(n1+n2);
//'
//' h2 = 0.001;
//' y1 <- y[1:n1]
//' X1 <- X[1:n1,]
//' y2 <- y0[(n1+1):(n1+n2)]
//' X2 <- X[(n1+1):(n1+n2),]
//' alpha0 <- 3
//' alpha <- 0.3
//' sz2 <- var(y2*alpha) * ((1-h2)/h2)
//' z <- alpha0 + y2*alpha + rnorm(n2,0,sqrt(sz2))
//' y = y1;
//'
//' mean.x1 = apply(X1,2,mean);
//' x1m = sweep(X1,2,mean.x1);
//' std.x1 = apply(x1m,2,sd)
//' x1p = sweep(x1m,2,std.x1,"/");
//' x1p = x1p/sqrt(dim(x1p)[2])
//'
//'  mean.x2 = apply(X2,2,mean);
//' x2m = sweep(X2,2,mean.x2);
//' std.x2 = apply(x2m,2,sd)
//' x2p = sweep(x2m,2,std.x2,"/");
//' x2p = x2p/sqrt(dim(x2p)[2])
//'
//' w2 = matrix(rep(1,n2),ncol=1);
//' w1 = matrix(rep(1,n1),ncol=1);
//' fm = CoMM_covar_pxem(y,z,x1p,x2p,w1,w2);
//'
//' @details
//' \code{CoMM_covar_pxem} fits the CoMM model using raw data.
//' @export
// [[Rcpp::export]]
List CoMM_covar_pxem(const arma::vec& y, const arma::vec& z, const arma::mat& x1,  const arma::mat& x2,  const arma::mat& w1,  const arma::mat& w2, const int constr = 0, const double epsStopLogLik = 1e-5, const int maxIter = 1000, const int pxem_indicator = 1){
    int n1 = y.n_elem, n2 = z.n_elem, p1 = x1.n_cols, p2 = x2.n_cols;
    if (p1 != p2){
        perror("The dimensions of x1 and x2 are not matched");
    }
    int p = p1;

    mat x1tx1 = x1.t()*x1, x2tx2 = x2.t()*x2, w1tw1 = w1.t()*w1, w2tw2 = w2.t()*w2, x1tw1 = x1.t()*w1, x2tw2 = x2.t()*w2;
    vec x1ty = x1.t()*y, x2tz = x2.t()*z, w2tz = w2.t()*z, w1ty = w1.t()*y;

    // initialization using lmm_pxem
    double sigma2y, sigma2beta, loglik0;
    vec beta0 =zeros<vec>(w1.n_cols);
    int iter0;
    mat Sigb = zeros<mat>(x1.n_cols,x1.n_cols);
    vec mub  = zeros<vec>(x1.n_cols);

    lmm_pxem_ptr2(y, w1, x1, maxIter,sigma2y,sigma2beta,beta0,loglik0,iter0,Sigb,mub);

    //declaration of variables used within loop
    mat Sigma_inv(p,p), R(p,p), invR(p,p), Sigma(p,p);
    vec mutmp(p), mu(p), res_y(n1), res_z(n2), x1mu(n1), x2mu(n2);
    double ztmp;
    mat I = eye<mat>(p,p);

    //initialization of parameters
    vec alpha0 = solve(w2tw2, w2tz);
    double alpha = 0, sigma2z = var(z), gam = 1, gam2 = gam*gam, alpha2 = alpha*alpha;

    vec loglik(maxIter), loglik2(maxIter);

    // initialization of likelihood
    loglik(0) = NAN;
    int Iteration = 1;
    for (int iter = 2; iter <= maxIter; iter ++ ) {
        // E-step
        Sigma_inv = x1tx1/sigma2y + x2tx2*(alpha2/sigma2z) + I /sigma2beta;
        //cout << "break : " << Sigma_inv << endl;
        R = chol(Sigma_inv);
        invR = inv(R);
        Sigma = invR*invR.t();
        mutmp = (x1ty - x1tw1*beta0)/sigma2y + (x2tz - x2tw2*alpha0)*(alpha/sigma2z);
        mu = Sigma*mutmp;

        //evaluate incomplete log-likelihood
        x1mu = x1*mu;
        x2mu = x2*mu;
        res_y = y - x1mu*gam - w1*beta0;
        res_z = z - alpha*x2mu - w2*alpha0;

        //loglik(iter - 1) = loglike_twas(R, res_y, res_z, mu, sigma2beta, sigma2y, sigma2z,n1,n2,p);
        loglike_twas(R, res_y, res_z, mu, sigma2beta, sigma2y, sigma2z,n1,n2,p, loglik(iter - 1));

        if ( loglik(iter - 1) - loglik(iter - 2) < -1e-7){
            perror("The likelihood failed to increase!");
        }
        // M-step
        mat trtmp1 = x1tx1 * Sigma;
        double tr1 = trace(trtmp1);

        if (pxem_indicator == 1){
            gam = sum( (x1ty - x1tw1*beta0) % mu )/(sum(x1mu % x1mu) + tr1);
            gam2 = gam*gam;
        }
        else {
            gam = 1; gam2 = 1;
        }

        beta0 = solve(w1tw1, w1ty - gam*x1tw1.t()*mu);
        alpha0= solve(w2tw2, w2tz - alpha*(x2tw2.t()*mu));

        ztmp = sum( (z - w2*alpha0)% x2mu);
        mat trtmp2 = x2tx2 * Sigma;
        double tr2 = trace(trtmp2);
        if (constr == 1){
            alpha = 0;
        }
        else {
            alpha = ztmp/(sum(x2mu % x2mu) + tr2);
        }
        alpha2 = alpha*alpha;

        res_y = y - x1mu*gam - w1*beta0;
        res_z = z - alpha*x2mu - w2*alpha0;

        sigma2y = (sum(res_y % res_y) + gam2 * tr1)/n1;
        sigma2z = (sum(res_z % res_z) + alpha2 * tr2)/n2;
        sigma2beta = (sum(mu % mu) + trace(Sigma))/p;

        // Reduction-step
        sigma2beta = gam2 * sigma2beta;
        alpha = alpha / gam;
        alpha2 = alpha*alpha;
        gam = 1;
        gam2 = 1;

        Iteration = iter;
        if (iter > 2){
            if (abs(loglik(iter - 1) - loglik(iter - 2)) < epsStopLogLik) {

                break;
            }
        }
    }


    vec loglik_out;
    int to = Iteration -1;
    loglik_out = loglik.subvec(0, to);

    double loglik_max = loglik(to);

    List output = List::create(Rcpp::Named("alpha") = alpha,
                               Rcpp::Named("alpha0") = alpha0,
                               Rcpp::Named("beta0") = beta0,
                               Rcpp::Named("sigma2y") = sigma2y,
                               Rcpp::Named("sigma2z") = sigma2z,
                               Rcpp::Named("sigma2beta") = sigma2beta,
                               Rcpp::Named("gam") = gam,
                               Rcpp::Named("loglik_seq") = loglik_out,
                               Rcpp::Named("loglik") = loglik_max,
                               Rcpp::Named("iteration") = Iteration-1);

    return output;
}
