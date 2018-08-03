#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

void lmm_pxem_jin(const arma::vec& y, const arma::mat& W, const arma::mat& X, const double& tol, const int& maxIter,
                  double& sigma2y, double& sigma2beta, arma::vec& beta0, double& loglik_max,
                  int& iteration, arma::mat& Sigb, arma::vec& mub){
  if (y.n_elem != X.n_rows || X.n_rows != W.n_rows){
    perror("The dimensions in outcome and covariates (X and W) are not matched");
  }

  int n = y.n_elem, p = X.n_cols;
  sigma2y = var(y), sigma2beta = sigma2y/p;

  if (beta0.n_elem != W.n_cols){
    perror("The dimensions in covariates are not matched in W and beta0");
  }

  if (p != mub.n_elem){
    perror("The dimensions in covariates are not matched in mub");
  }

  if (p != Sigb.n_cols){
    perror("The dimensions in covariates are not matched in Sigb");
  }

  mat XtX = X.t()*X, WtW = W.t()*W, XtW = X.t()*W;
  vec Xty = X.t()*y, Wty = W.t()*y;
  double gam, gam2;  // parameter expansion

  vec dd;
  mat uu;

  eig_sym(dd, uu, XtX);

  mat wtxu = W.t()*X*uu; //q-by-p

  vec dd2;
  mat uu2;

  mat XXt = X*X.t();

  eig_sym(dd2, uu2, XXt);

  vec Cm(p), Cs(p), Cm2(p);
  vec loglik(maxIter);
  vec utxty(p), utxty2(p), utxty2Cm(p), utxty2Cm2(p);
  vec yt(n), yt2(n);

  // evaluate loglik at initial values
  vec uty = uu2.t()*(y - W * beta0);

  vec tmpy = uty / sqrt(dd2 * sigma2beta + sigma2y);
  vec tmpy2 = pow(tmpy,2);
  loglik(0) = -0.5 * sum(log(dd2 * sigma2beta + sigma2y)) - 0.5 * sum(tmpy2);

  int Iteration = 1;
  for (int iter = 2; iter <= maxIter; iter ++ ) {
    // E-step
    Cm = sigma2y / sigma2beta +  dd;
    Cs = 1 / sigma2beta +  dd / sigma2y;
    Cm2 = pow(Cm , 2);
    // M-step
    utxty = uu.t() * (Xty - XtW * beta0); // p-by-1
    utxty2 = pow(utxty, 2);

    utxty2Cm = utxty % utxty / Cm;
    utxty2Cm2 = utxty2Cm/Cm;

    gam = sum(utxty2Cm) / ( sum(dd % utxty2Cm2) + sum(dd / Cs));
    gam2 = pow(gam , 2);

    //sigma2beta = ( sum(utxty.t() * diagmat(1 / Cm2) * utxty) + sum(1 / Cs)) / p;
    sigma2beta = ( sum(utxty2Cm2) + sum(1 / Cs)) / p;

    yt = y - W*beta0;
    yt2 = pow(yt , 2);
    sigma2y = (sum(yt2) - 2 * sum(utxty2Cm) * gam + sum(utxty2Cm2 % dd) * gam2 + gam2 * sum(dd / Cs)) / n;

    beta0 = solve(WtW, Wty - gam2 * (wtxu * (utxty / Cm)));

    //reduction and reset
    sigma2beta = gam2 * sigma2beta;

    //evaluate loglik
    uty = uu2.t()*(y - W * beta0);
    tmpy = uty / sqrt(dd2 * sigma2beta + sigma2y);
    tmpy2 = pow(tmpy,2);
    loglik(iter - 1) = - 0.5 * sum(log(dd2 * sigma2beta + sigma2y)) - 0.5 * sum(tmpy2);

    if ( loglik(iter - 1) - loglik(iter - 2) < 0 ){
      perror("The likelihood failed to increase!");
    }

    Iteration = iter;
    if (abs(loglik(iter - 1) - loglik(iter - 2)) < 1e-5) {

      break;
    }
  }

  Cm = sigma2y / sigma2beta + dd;
  Cs = 1 / sigma2beta + dd / sigma2y;
  Sigb = uu*diagmat(1 / Cs)*uu.t();
  mub = uu * (utxty / Cm);

  vec loglik_out;
  int to = Iteration -1;
  loglik_out = loglik.subvec(0, to);

  loglik_max = loglik(to);
  iteration = Iteration -1;
}

void lmm_pxem_ptr2(const arma::vec& y, const arma::mat& W, const arma::mat& X,  const int& maxIter,
              double& sigma2y, double& sigma2beta, arma::vec& beta0, double& loglik_max,
              int& iteration, arma::mat& Sigb, arma::vec& mub){

  int n = y.n_elem, p = X.n_cols;

  if (y.n_elem != X.n_rows || X.n_rows != W.n_rows){
    perror("The dimensions in outcome and covariates (X and W) are not matched");
  }

  if (beta0.n_elem != W.n_cols){
    perror("The dimensions in covariates are not matched in W and beta0");
  }

  if (p != mub.n_elem){
    perror("The dimensions in covariates are not matched in mub");
  }

  if (p != Sigb.n_cols){
    perror("The dimensions in covariates are not matched in Sigb");
  }

  mat XtX = X.t()*X, WtW = W.t()*W, WtX = W.t()*X;
  vec Xty = X.t()*y, Wty = W.t()*y;

  vec SWy;
  mat SWX;

  if(W.n_cols==1){
    SWy = mean(y);
    SWX = mean(X,0);
  } else{
    SWy = solve(WtW, Wty);
    SWX = solve(WtW, WtX);
  }


  double gam, gam2;  // parameter expansion

  vec eVal;
  mat eVec;

  eig_sym(eVal, eVec, XtX);

  // initialize
  sigma2y = var(y);
  sigma2beta = sigma2y/p;
  beta0 = SWy - SWX * mub;
  vec loglik(maxIter);
  loglik(0) = -datum::inf;

  vec D;
  vec Xmu;
  vec y_bar = y - W * beta0;
  double y_Xmu2, E;

  iteration = maxIter-1;
  for (int iter = 1; iter < maxIter; iter ++ ) {
    // E-step
    D = 1 / sigma2beta +  eVal / sigma2y;
    mub = 1/sigma2y * eVec * (eVec.t() * (X.t() * y_bar) / D);
    Xmu = X * mub;
    y_Xmu2 = sum(pow(y_bar-Xmu,2));

    // Evaluate loglik
    E = y_Xmu2/(2*sigma2y) + accu(pow(mub,2))/(2*sigma2beta);
    loglik(iter) = - p*log(sigma2beta)/2 - n*log(sigma2y)/2 - E - sum(log(D))/2 - n/2*log(2*datum::pi);

    if ( loglik(iter) - loglik(iter - 1) < 0 ){
      perror("The likelihood failed to increase!");
    }

    if (abs(loglik(iter) - loglik(iter - 1)) < 1e-10) {
      iteration = iter;
      break;
    }

    // M-step
    gam = sum(y_bar % Xmu) / (accu(pow(Xmu,2)) + sum(eVal/D));
    gam2 = pow(gam , 2);

    beta0 = SWy - (SWX * mub) * gam;
    y_bar = y - W * beta0;;

    sigma2y = sum(pow(y_bar-Xmu*gam,2))/n + gam2 * sum(eVal/D)/n;
    sigma2beta = accu(pow(mub,2))/p + sum(1/D)/p;

    // Reduction step
    sigma2beta = gam2 * sigma2beta;
    // gam = 1;
    // gam2 = pow(gam , 2);
  }

  vec loglik_out;
  loglik_out = loglik.subvec(0, iteration);

  loglik_max = loglik(iteration);
}

// [[Rcpp::export]]
Rcpp::List lmm_pxem2(const arma::vec y, const arma::mat w, const arma::mat x, const int maxIter){

    double sigma2y, sigma2beta, loglik;
    vec beta0 =zeros<vec>(w.n_cols);
    int iter;
    mat Sigb = zeros<mat>(x.n_cols,x.n_cols);
    vec mub  = zeros<vec>(x.n_cols);

    lmm_pxem_ptr2(y, w, x, maxIter,sigma2y,sigma2beta,beta0,loglik,iter,Sigb,mub);

    List output = List::create(Rcpp::Named("sigma2y") = sigma2y,
                               Rcpp::Named("sigma2beta") = sigma2beta,
                               Rcpp::Named("beta0") = beta0,
                               //Rcpp::Named("loglik_seq") = loglik_out,
                               Rcpp::Named("loglik") = loglik,
                               Rcpp::Named("iteration") = iter,
                               Rcpp::Named("Sigb") = Sigb,
                               Rcpp::Named("mub") = mub);

    return output;

}
