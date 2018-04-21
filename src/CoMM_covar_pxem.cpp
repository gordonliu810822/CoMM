#include <RcppArmadillo.h>
//#include <armadillo>
//#include <Rcpp.h>
//#include <omp.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

//const double PI  = 3.141592653589793238463;

double loglike_twas(arma::mat R, arma::vec res_y, arma::vec res_z, arma::vec mu, double sigma2beta, double sigma2y, double sigma2z, int n1,
					int n2, int p){
	double Eab, loglik;
	Eab = sum(res_y % res_y)/2/sigma2y + sum(res_z % res_z)/2/sigma2z + sum(mu % mu)/2/sigma2beta;
	/*cout << "Eab: " << Eab << endl;
	cout << "sum(log(R.diag())): " << sum(log(R.diag())) << endl;
	cout << "p /2 *log(sigma2beta): " << log(sigma2beta) << endl;
	cout << "n1/2*log(sigma2y): " << n1 / 2 * log(sigma2y) << endl;
	cout << "n2/2*log(sigma2z): " << n2 / 2 * log(sigma2z) << endl;*/

	loglik = -log(sigma2beta)*p*0.5 - log(sigma2y)*n1*0.5 - log(sigma2z)*n2 *0.5 - sum(log(R.diag())) - Eab; //- (n1+n2+p)/2*log(2*M_PI)

	return loglik;
}

double loglike_twas2(mat w1, mat w2, vec alpha0, vec beta0, mat x1x1t, mat x1x2t, mat x2x2t,
					 double alpha, double sigma2beta, double sigma2y, double sigma2z, vec y, vec z, int n1, 
					int n2, int p){
	vec hmean(n1 + n2);
	hmean.subvec(0,n1-1) = w1*beta0;
	hmean.subvec(n1,n1+n2-1) = w2*alpha0;
	mat Sigma11 = x1x1t*sigma2beta + sigma2y*eye<mat>(n1,n1), Sigma12 = x1x2t*sigma2beta*alpha;
	mat Sigma22 = x2x2t*sigma2beta*alpha*alpha + sigma2z*eye<mat>(n2,n2);
	mat tmp1 = join_rows(Sigma11, Sigma12), tmp2 = join_rows(Sigma12.t(), Sigma22) ;
	mat Sigma = join_cols(tmp1,tmp2);
	mat R = chol(Sigma), invR = inv(R);
	double logdetSigma = sum(log(R.diag()));

	vec zz = join_cols(y,z);
	vec tmp = invR.t()*(zz-hmean);
	double loglik = - logdetSigma - 0.5* sum(tmp % tmp);

	return loglik;
}

//' @title
//' CoMM
//' @description
//' fit CoMM for a single gene
//'
//' @param y  a vector for the expression of a gene.
//' @param z  a vector for the phenotype of GWAS.
//' @param X1  a standardized genotype matrix for eQTL data.
//' @param X2  a standardized genotype matrix for GWAS data.
//' @param w1  a matrix of coveriates for the expression of a gene.
//' @param w2  a matrix of coveriates for the GWAS.
//' @param sigma2beta  a initial value of sigma2beta.
//' @param sigma2y  a initial value of sigma2y.
//' @param beta0  a initial value of beta0.
//' @param constr  indicator value for constraint. When constr = 1, alpha = 0, otherwise, alpha is free to vary.
//' @param epsStopLogLik convergence criteria. default is 1e-5.
//' @param maxIter maximum iteraion number. default is 1000.
//'
//' @return List of model parameters
//'
//' @export
// [[Rcpp::export]]
List CoMM_covar_pxem(arma::vec y, arma::vec z, arma::mat x1,  arma::mat x2,  arma::mat w1,  arma::mat w2, double sigma2beta,
					 double sigma2y, arma::vec beta0, int constr, double epsStopLogLik=1e-5, int maxIter=1000){
	int n1 = y.n_elem, n2 = z.n_elem, p1 = x1.n_cols, p2 = x2.n_cols; //q1 = w1.n_cols, q2 = w2.n_cols;
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
	vec alpha0 = solve(w2tw2, w2tz);
	double alpha = 0, sigma2z = var(z), gam = 1, gam2 = gam*gam, alpha2 = alpha*alpha;

	//mat x1x1t = x1*x1.t(), x1x2t = x1*x2.t(), x2x2t = x2*x2.t();

	vec loglik(maxIter), loglik2(maxIter);

	/*// initialization of likelihood
	res_y = y - x1mu*gam - w1*beta0;
	res_z = z - alpha*x2mu - w2*alpha0;
	Sigma_inv = x1tx1 / sigma2y + x2tx2*(alpha2 / sigma2z) + I / sigma2beta;
	R = chol(Sigma_inv);
	invR = inv(R);
	Sigma = invR*invR.t();
	mutmp = (x1ty - x1tw1*beta0) / sigma2y + (x2tz - x2tw2*alpha0)*(alpha / sigma2z);
	mu = Sigma*mutmp;

	loglik(0) = loglike_twas(R, res_y, res_z, mu, sigma2beta, sigma2y, sigma2z, n1, n2, p);*/
	loglik(0) = NAN;
	//printf("1-st iter: loglik = %11.8f \n", loglik(0));

	/*loglik2(0)= loglike_twas2(w1, w2, alpha0, beta0, x1x1t, x1x2t, x2x2t,
					  alpha, sigma2beta, sigma2y, sigma2z, y, z, n1, n2, p);*/
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

		loglik(iter - 1) = loglike_twas(R, res_y, res_z, mu, sigma2beta, sigma2y, sigma2z,n1,n2,p);
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

	//cout << "Finished in " << Iteration << ", loglik = " << loglik(Iteration -1) << endl;

	vec loglik_out;
	int to = Iteration -1;
	loglik_out = loglik.subvec(0, to);

	List output = List::create(Rcpp::Named("alpha") = alpha,
		Rcpp::Named("alpha0") = alpha0,
		Rcpp::Named("beta0") = beta0,
		Rcpp::Named("sigma2y") = sigma2y,
		Rcpp::Named("sigma2z") = sigma2z,
		Rcpp::Named("sigma2beta") = sigma2beta,
		Rcpp::Named("R") = R,
		Rcpp::Named("mu") = mu,
		Rcpp::Named("gam") = gam,
		Rcpp::Named("loglik") = loglik_out);//.subvec(0,Iteration-1));

	return output;
}

