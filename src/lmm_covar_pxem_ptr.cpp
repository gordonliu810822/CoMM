#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

void lmm_pxem_ptr(const arma::vec& y, const arma::mat& w, const arma::mat& x, const int& maxIter,
			  double& sigma2y, double& sigma2beta, arma::vec& beta0, double& loglik_max,
			  int& iteration, arma::mat& Sigb, arma::vec& mub){
	if (y.n_elem != x.n_rows || x.n_rows != w.n_rows){
		perror("The dimensions in outcome and covariates (x and w) are not matched");
	}

	sigma2y = var(y), sigma2beta = sigma2y;
    int n = y.n_elem, p = x.n_cols;
	if (beta0.n_elem != w.n_cols){
		perror("The dimensions in covariates are not matched in w and beta0");
	}

	if (p != mub.n_elem){
		perror("The dimensions in covariates are not matched in mub");
	}

	if (p != Sigb.n_cols){
		perror("The dimensions in covariates are not matched in Sigb");
	}

	mat xtx = x.t()*x, wtw = w.t()*w, xtw = x.t()*w;
	vec xty = x.t()*y, wty = w.t()*y;
	double gam, gam2;  // parameter expansion

	vec  dd;
	mat uu;

	eig_sym(dd, uu, xtx);

	mat wtxu = w.t()*x*uu; //q-by-p

	vec dd2;
	mat uu2;

	mat tmp = x*x.t();

	eig_sym(dd2, uu2, tmp);

	vec Cm(p), Cs(p), Cm2(p);
	vec loglik(maxIter);
	vec utxty(p), utxty2(p), utxty2Cm(p), utxty2Cm2(p);
	vec yt(n), yt2(n);

	// evaluate loglik at initial values
	vec uty = uu2.t()*(y - w * beta0);

	vec tmpy = uty / sqrt(dd2 * sigma2beta + sigma2y);
	vec tmpy2 = pow(tmpy,2);
	loglik(0) = -0.5 * sum(log(dd2 * sigma2beta + sigma2y)) - 0.5 * sum(tmpy2);

	//cout << "initial :sigma2beta: " << sigma2beta <<"; sigma2y: " << sigma2y << "; loglik: " << loglik(0) << endl;
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
		//cout << " iter: " << iter << "; 1/Cs: " << sum(1 / Cs) << ";utxty2Cm2: " << sum(utxty2Cm2) <<"; gam: " << sum(gam) <<  "; utxty2Cm: " <<	sum(utxty2Cm) <<endl;

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

		//cout << "sigma2beta: " << sigma2beta <<"; sigma2y: " << sigma2y << "; beta0: " << beta0 << "; loglik: " <<	loglik(iter - 1) << "; sum(tmpy2): " << sum(tmpy2) <<endl;
		if ( loglik(iter - 1) - loglik(iter - 2) < 0 ){
			perror("The likelihood failed to increase!");
		}

		Iteration = iter;
		if (abs(loglik(iter - 1) - loglik(iter - 2)) < 1e-10) {

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

void lmm_pxem_test(const arma::vec& y, const arma::mat& w, const arma::mat& x, const int& maxIter, double& sigma2y, double& sigma2beta, arma::vec& beta0, double& loglik_max, int& iteration, arma::vec& mub){
	// using eigen decomposion on xxt (may be semi-positive-definite): finally fixed!
    if (y.n_elem != x.n_rows || x.n_rows != w.n_rows){
        perror("The dimensions in outcome and covariates (x and w) are not matched");
    }

    sigma2y = var(y), sigma2beta = sigma2y;
    int n = y.n_elem, p = x.n_cols;
    if (beta0.n_elem != w.n_cols){
        perror("The dimensions in covariates are not matched in w and beta0");
    }

    if (p != mub.n_elem){
        perror("The dimensions in covariates are not matched in mub");
    }

    //if (p != Sigb.n_cols){
    //    perror("The dimensions in covariates are not matched in Sigb");
    //}

    mat wtw = w.t()*w, xtw = x.t()*w;
    vec xty = x.t()*y, wty = w.t()*y;
    double gam, gam2;  // parameter expansion

    //vec  dd;
    //mat uu;

    //eig_sym(dd, uu, xtx);

    vec dd3; //n-by-1
    mat uu3;

    mat tmp = x*x.t();

    eig_sym(dd3, uu3, tmp);

    //mat uu(p,n);

	mat uu2 = uu3.cols(find(dd3>1e-5)); //n-by-n0
	vec dd2 = dd3.elem(find(dd3>1e-5));//n0-by-1

	uword n0 = dd2.n_elem;

	mat D = repmat(sqrt(dd2.t()),p,1);

    mat uu = x.t()*uu2/D;  //p-by-n0, uu*dd2*uut = xtx

    mat wtxu = w.t()*x*uu; //q-by-n0

    //cout << "n: " << n << "; n0: " << n0 << "; sum(dd): " << sum(dd2) << endl;

    vec Cm(n0), Cs(n0);// Cm2(p);
    vec loglik(maxIter);
    vec utxty(n0), utxty2(n0), utxty2Cm(n0), utxty2Cm2(n0);
    vec yt(n0), yt2(n0);
    //cout << "break 1..." << endl;

    // evaluate loglik at initial values
    vec uty = uu2.t()*(y - w * beta0); //n0-by-1

    vec tmpy = uty / sqrt(dd2 * sigma2beta + sigma2y); //n0-by-1
    vec tmpy2 = pow(tmpy,2);  //n0-by-1
    loglik(0) = -0.5 * sum(log(dd2 * sigma2beta + sigma2y)) - 0.5 * sum(tmpy2);
	//cout << "initial :sigma2beta: " << sigma2beta <<"; sigma2y: " << sigma2y << "; loglik: " << loglik(0) << endl;

    //if (constr == 0){
        int Iteration = 1;
        for (int iter = 2; iter <= maxIter; iter ++ ) {
            // E-step
            Cm = sigma2y / sigma2beta +  dd2;
            Cs = 1 / sigma2beta +  dd2 / sigma2y;
            //Cm2 = pow(Cm , 2);
            // M-step
            utxty = uu.t() * (xty - xtw * beta0); // n0-by-1
            utxty2 = pow(utxty, 2);// n0-by-1

            utxty2Cm = utxty % utxty / Cm;// n0-by-1
            utxty2Cm2 = utxty2Cm/Cm;// n0-by-1

            gam = sum(utxty2Cm) / ( sum(dd2 % utxty2Cm2) + sum(dd2 / Cs));
            gam2 = pow(gam , 2);
            //cout << " iter: " << iter << "; 1/Cs: " << sum(1 / Cs) << ";utxty2Cm2: " << sum(utxty2Cm2) <<"; gam: " << sum(gam) <<  "; utxty2Cm: " <<    sum(utxty2Cm) <<endl;

            //sigma2beta = ( sum(utxty.t() * diagmat(1 / Cm2) * utxty) + sum(1 / Cs)) / p;
            sigma2beta = ( sum(utxty2Cm2) + sum(1 / Cs) + sigma2beta*(p-n0)) / p;

            yt = y - w*beta0;
            yt2 = pow(yt , 2);
            sigma2y = (sum(yt2) - 2 * sum(utxty2Cm) * gam + sum(utxty2Cm2 % dd2) * gam2 + gam2 * sum(dd2 / Cs)) / n;

            beta0 = solve(wtw, wty - gam2 * (wtxu * (utxty / Cm)));

            //reduction and reset
            sigma2beta = gam2 * sigma2beta;

            //evaluate loglik
            uty = uu2.t()*(y - w * beta0);
            tmpy = uty / sqrt(dd2 * sigma2beta + sigma2y);
            tmpy2 = pow(tmpy,2);
            loglik(iter - 1) = - 0.5 * sum(log(dd2 * sigma2beta + sigma2y) ) - 0.5*log(sigma2y)*(n-n0) - 0.5 * sum(tmpy2);
            //cout << "sigma2beta: " << sigma2beta <<"; sigma2y: " << sigma2y << "; beta0: " << beta0 << "; loglik: " <<    loglik(iter - 1) << "; sum(tmpy2): " << sum(tmpy2) <<endl;

            if ( loglik(iter - 1) - loglik(iter - 2) < -1e-10 ){
                perror("The likelihood failed to increase!");
            }

            Iteration = iter;
            if (abs(loglik(iter - 1) - loglik(iter - 2)) < 1e-10) {

                break;
            }
        }
    

        Cm = sigma2y / sigma2beta + dd2;
        Cs = 1 / sigma2beta + dd2 / sigma2y;
        //Sigb = uu*diagmat(1 / Cs)*uu.t();
        mub = uu * (utxty / Cm);

        vec loglik_out;
        int to = Iteration -1;
        loglik_out = loglik.subvec(0, to);

        loglik_max = loglik(to);
        iteration = Iteration -1;
    //}
    
    /*if (constr == 1){
        beta0 = solve(wtw, wty);
        vec res = y - w * beta0;
        sigma2y = sum(res % res)/n;
        // sigma2beta = 0;
        uty = uu2.t()*res;
        
        vec tmpy = uty / sqrt(sigma2y); //n0-by-1
        vec tmpy2 = pow(tmpy,2);  //n0-by-1
        loglik_max = -0.5 * sum(log(sigma2y)) - 0.5*log(sigma2y)*(n-n0) - 0.5 * sum(tmpy2);
    }*/
}

//' @author Jin Liu, \email{jin.liu@duke-nus.edu.sg}
//' @title
//' CoMM
//' @description
//' CoMM to dissecting genetic contributions to complex traits by leveraging regulatory information.
//'
//' @param y  gene expression vector.
//' @param x  normalized genotype (cis-SNPs) matrix for eQTL.
//' @param w  covariates file for eQTL data.
//' @param maxIter  maximum iteration (default is 1000).
//'
//' @return List of model parameters
//'
//' @examples
//' L = 1; M = 100; rho =0.5
//' n1 = 350; n2 = 5000;
//' X <- matrix(rnorm((n1+n2)*M),nrow=n1+n2,ncol=M);
//'
//' beta_prop = 0.2;
//' b = numeric(M);
//' m = M * beta_prop;
//' b[sample(M,m)] = rnorm(m);
//' h2y = 0.05;
//' b0 = 6;
//'
//' y0 <- X%*%b + b0;
//' y  <- y0 + (as.vector(var(y0)*(1-h2y)/h2y))^0.5*rnorm(n1+n2);
//' h2 = 0.001;
//' y1 <- y[1:n1]
//' X1 <- X[1:n1,]
//' y = y1;
//'
//' mean.x1 = apply(X1,2,mean);
//' x1m = sweep(X1,2,mean.x1);
//' std.x1 = apply(x1m,2,sd)
//' x1p = sweep(x1m,2,std.x1,"/");
//' x1p = x1p/sqrt(dim(x1p)[2])
//' w1 = matrix(rep(1,n1),ncol=1);
//'
//' fm0 = lmm_pxem(y, w1,x1p, 100)
//'
//' @details
//' \code{lmm_pxem} fits the linear mixed model. (n < p)
//' @export
// [[Rcpp::export]]
Rcpp::List lmm_pxem(const arma::vec y, const arma::mat w, const arma::mat x, const int maxIter){

    double sigma2y, sigma2beta, loglik;
    vec beta0 =zeros<vec>(w.n_cols);
    int iter;
    mat Sigb = zeros<mat>(x.n_cols,x.n_cols);
    vec mub  = zeros<vec>(x.n_cols);
    //mat uu;// = zeros<mat>(x.n_cols,x.n_cols);
	//mat D;// = zeros<mat>(x.n_cols,x.n_rows);
    lmm_pxem_ptr(y, w, x, maxIter,sigma2y,sigma2beta,beta0,loglik,iter,Sigb,mub);

    //lmm_pxem_test(y, w, x, maxIter,sigma2y,sigma2beta,beta0,loglik,iter,mub);

    List output = List::create(Rcpp::Named("sigma2y") = sigma2y,
                               Rcpp::Named("sigma2beta") = sigma2beta,
                               Rcpp::Named("beta0") = beta0,
                               //Rcpp::Named("loglik_seq") = loglik_out,
                               Rcpp::Named("loglik") = loglik,
                               Rcpp::Named("iteration") = iter,
                               Rcpp::Named("Sigb") = Sigb,
                               Rcpp::Named("mub") = mub);
                               //Rcpp::Named("uu") = uu,
							   //Rcpp::Named("D") = D);

    return output;

}

// [[Rcpp::export]]
Rcpp::List lmm_pxem_test(const arma::vec y, const arma::mat w, const arma::mat x, const int maxIter){
    // constr= 0, 1; 0 for sigma2beta !=0 and 1 for sigma2beta = 0;
    double sigma2y, sigma2beta, loglik;
    vec beta0 =zeros<vec>(w.n_cols);
    int iter;
    //mat Sigb = zeros<mat>(x.n_cols,x.n_cols);
    vec mub  = zeros<vec>(x.n_cols);
    //mat uu;// = zeros<mat>(x.n_cols,x.n_cols);
	//mat D;// = zeros<mat>(x.n_cols,x.n_rows);

    lmm_pxem_test(y, w, x, maxIter,sigma2y,sigma2beta,beta0,loglik,iter,mub);

    List output = List::create(Rcpp::Named("sigma2y") = sigma2y,
                               Rcpp::Named("sigma2beta") = sigma2beta,
                               Rcpp::Named("beta0") = beta0,
                               //Rcpp::Named("loglik_seq") = loglik_out,
                               Rcpp::Named("loglik") = loglik,
                               Rcpp::Named("iteration") = iter,
                               //Rcpp::Named("Sigb") = Sigb,
                               Rcpp::Named("mub") = mub);
                               //Rcpp::Named("uu") = uu,
							   //Rcpp::Named("D") = D);

    return output;

}
