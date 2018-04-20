#include <RcppArmadillo.h>
//#include <Rcpp.h>
//#include <omp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
List lmm_pxem(arma::vec y, arma::mat w, arma::mat x, int maxIter){
	// y: outcome; w: design matrix of fixed effects; x: matrix of random effects;
	if (y.n_elem != x.n_rows || x.n_rows != w.n_rows){
		perror("The dimensions in outcome and covariates (x and w) are not matched");
	}

	double sigma2y = var(y), sigma2beta = sigma2y;
	int n = y.n_elem, p = x.n_cols, q = w.n_cols;
	//cout << "break :sigma2y = " << sigma2y << "; sigma2beta = " << sigma2beta << endl;

	vec beta0 = zeros<vec>(q);
	mat xtx = x.t()*x, wtw = w.t()*w, xtw = x.t()*w;
	vec xty = x.t()*y, wty = w.t()*y;
	double gam, gam2;  // parameter expansion

	vec  dd;
	mat uu;

	//cout << "break :col sum of x :" << sum( x, 0 ) << endl;
	//cout << "break :x : " << x.n_rows << "-by-" << x.n_cols << endl;
	//cout << "break :w : " << w.n_rows << "-by-" << w.n_cols << endl;
	//cout << "break :xtx : " << xtx.n_rows << "-by-" << xtx.n_cols << endl;
	//cout << "break :xtx : " << xtx << endl;
	eig_sym(dd, uu, xtx);

	//cout << "break :col sum of uu " << sum( uu, 0 ) << endl;
	mat wtxu = w.t()*x*uu; //q-by-p
	//cout << "break :wtxu : " << wtxu.n_rows << "-by-" << wtxu.n_cols << endl;

	vec dd2;
	mat uu2;

	mat tmp = x*x.t();

	//cout << "break :xxt : " << tmp.n_rows << "-by-" << tmp.n_cols << endl;
	eig_sym(dd2, uu2, tmp);

	vec Cm(p), Cs(p), Cm2(p);
	vec loglik(maxIter);
	vec utxty(p), utxty2(p), utxty2Cm(p), utxty2Cm2(p);
	vec yt(n), yt2(n);

	// evaluate loglik at initial values
	vec uty = uu2.t()*(y - w * beta0);

	vec tmpy = uty / sqrt(dd2 * sigma2beta + sigma2y);
	vec tmpy2 = pow(tmpy,2);
	//cout << "break : tmpy2: " << - 0.5 * sum(tmpy2) << endl;
	//cout << "break : dd2: " << sum(dd2) << endl;
	//cout << "break : tmp: " << -0.5 * sum(log(dd2 % sigma2beta_vec + sigma2y)) << endl;
	loglik(0) = -0.5 * sum(log(dd2 * sigma2beta + sigma2y)) - 0.5 * sum(tmpy2);
	//cout << "break : loglik(1)" << loglik(0) << endl;

	int Iteration = 1;
	for (int iter = 2; iter <= maxIter; iter ++ ) {
		// E-step
		Cm = sigma2y / sigma2beta +  dd;
		Cs = 1 / sigma2beta +  dd / sigma2y;
		Cm2 = pow(Cm , 2);
		// M-step
		utxty = uu.t() * (xty - xtw * beta0); // p-by-1
		utxty2 = pow(utxty, 2);

		utxty2Cm = utxty % utxty / Cm;
		utxty2Cm2 = utxty2Cm/Cm;

		gam = sum(utxty2Cm) / ( sum(dd % utxty2Cm2) + sum(dd / Cs));
		gam2 = pow(gam , 2);

		//sigma2beta = ( sum(utxty.t() * diagmat(1 / Cm2) * utxty) + sum(1 / Cs)) / p;
		sigma2beta = ( sum(utxty2Cm2) + sum(1 / Cs)) / p;

		yt = y - w*beta0;
		yt2 = pow(yt , 2);
		sigma2y = (sum(yt2) - 2 * sum(utxty2Cm) * gam + sum(utxty2Cm2 % dd) * gam2 + gam2 * sum(dd / Cs)) / n;

		beta0 = solve(wtw, wty - gam2 * (wtxu * (utxty / Cm)));
		
		//reduction and reset
		sigma2beta = gam2 * sigma2beta;

		//evaluate loglik
		uty = uu2.t()*(y - w * beta0);
		tmpy = uty / sqrt(dd2 * sigma2beta + sigma2y);
		tmpy2 = pow(tmpy,2);
		loglik(iter - 1) = - 0.5 * sum(log(dd2 * sigma2beta + sigma2y)) - 0.5 * sum(tmpy2);
		// cout << iter << "-th iter: loglik = " << loglik(iter - 1) << "; sigma2y "  << sigma2y << "; sigma2beta " << sigma2beta << "; diff : " << loglik(iter - 1) - loglik(iter - 2) << endl;
		// cout << iter << "-th iter: loglik = " << loglik(iter - 1) << "; diff : " << loglik(iter - 1) - loglik(iter - 2) << endl;

		if ( loglik(iter - 1) - loglik(iter - 2) < 0 ){
			perror("The likelihood failed to increase!");
		}

		Iteration = iter;
		if (abs(loglik(iter - 1) - loglik(iter - 2)) < 1e-6) {
			
			break;
		}
	}
	
	Cm = sigma2y / sigma2beta + dd;
	Cs = 1 / sigma2beta + dd / sigma2y;
	//utxty = uu.t() * (xty - xtw * beta0);

	mat Sigb = uu*diagmat(1 / Cs)*uu.t();
	
	vec mub1 = uu * (utxty / Cm);

	/*mat I = eye<mat>(p,p); 
	mat Sigma_inv = xtx / sigma2y + I / sigma2beta;
	//mat tmpinv = xtx + eye<mat>(p, p) *sigma2y / sigma2beta;
	//vec mub2 = solve(tmpinv, (xty - xtw * beta0));
	mat R = chol(Sigma_inv);
	mat invR = inv(R);
	mat Sigma = invR*invR.t();
	vec mu = (xty - xtw*beta0)/sigma2y;
	mu = Sigma*mu;*/
	// mub2 = invSigb2*(xty - xtw * beta0) / sigma2y;
	// mat Sigb = inv(invSigb2);

	vec loglik_out;
	int to = Iteration -1;
	loglik_out = loglik.subvec(0, to);

	List output = List::create(Rcpp::Named("sigma2y") = sigma2y,
		Rcpp::Named("sigma2beta") = sigma2beta,
		Rcpp::Named("beta0") = beta0,
		Rcpp::Named("loglik") = loglik_out,
		//Rcpp::Named("invSigb") = invSigb,
		//Rcpp::Named("invSigb") = Sigma_inv,
		//Rcpp::Named("Sigb") = Sigma,
		Rcpp::Named("Sigb") = Sigb,
		Rcpp::Named("mub") = mub1);
		//Rcpp::Named("mub2") = mub2);//.subvec(0,Iteration-1));

	return output;
}
