#include "RcppArmadillo.h"
#include "math.h"
#include <cassert>
#include <cmath>
#include "SCLMM.hpp"
#include "lmm_covar_pxem_ptr.hpp"

using namespace std;
using namespace arma;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R



double LLD(mat muinSRinSmu, mat xtx, vec hatmus2, mat inSRinS, vec mu, double alphag, double sigma2mu, int Mg, double sigma2e, double term1_1, vec v2){

	mat SIGMA = inv(xtx / sigma2e + alphag*alphag * inSRinS + 1 / sigma2mu*eye(Mg, Mg));
	double term1 = term1_1 - 0.5 / sigma2e*trace(SIGMA*xtx);
	double term2 = as_scalar(alphag*sum(hatmus2%mu) - 0.5*alphag*alphag*muinSRinSmu - 0.5*alphag*alphag*trace(SIGMA*inSRinS));
	double term3 = -0.5*Mg*sum(log(sigma2mu)) - 0.5*(sum(mu%mu) + trace(SIGMA)) / sigma2mu;
	double term4 = sum(log(diagvec(chol(SIGMA))));
	double out = term1 + term2 + term3 + term4;
	return out;
}


ObjSCLMM rcpparma_SCLMM_IS(mat &xr, vec &yr, mat &Wr, vec &hatmur, vec &hatsr, mat &R, Options_SCLMM* opts, bool px = 0){

	// check number of input arguments
	bool dispF = opts->dispF;
	int dispEvery = opts->display_gap;
	int max_iter = opts->max_iter;
	double tol = opts->epsStopLogLik;
	bool fix_alphag = opts->fix_alphag;

	double diff = datum::inf;
	double Fold = -datum::inf;

	int n = xr.n_rows;
	int Mg = hatmur.n_rows;
	int q = Wr.n_cols;

	mat x(xr.begin(), n, Mg, false);
	vec y(yr.begin(), yr.size(), false);
	mat W(Wr.begin(), n, q, false);
	vec hatmu(hatmur.begin(), hatmur.size(), false);
	vec hats(hatsr.begin(), hatsr.size(), false);

	double sigma2y, sigma2beta, loglik;
	vec beta0 = zeros<vec>(W.n_cols);
	int iterLmm, maxIter = 100000;
	mat Sigb = zeros<mat>(x.n_cols, x.n_cols);
	vec mub = zeros<vec>(x.n_cols);
	lmm_pxem_ptr(y, W, x, maxIter, sigma2y, sigma2beta, beta0, loglik, iterLmm, Sigb, mub);

	double gam = 1;
	bool px_ = px;
	int iter = 1;
	colvec mu = zeros(Mg, 1); // initial variational mean
	rowvec lowerbound = zeros(1, max_iter + 1);

	double term1, term2, term3, term4, term1_1;
	double F = 0;
	rowvec Pr = zeros(1, 4);

	mat SIGMA = zeros(Mg, Mg);
	mat S = diagmat(hats);
	vec v2 = zeros(Mg, 1);
	vec beta = zeros(q, 1);


	// predefined constant
	mat xtx = x.t()*x;
	vec diagxtx = sum(x%x).t();
	vec hatmus2 = hatmu / hats / hats;
	mat inS = inv_sympd(S);
	mat RinS = R*inS;
	mat inSRinS = eye(Mg, Mg)*inS*R*inS;
	vec diaginSRinS = diagvec(inSRinS);
	vec RinSmu = RinS*mu;
	vec ytilde = x*mu;
	vec ytilde_j;
	mat invSRS = inv(S*R*S);

	mat muinSRinSmu;
	vec xxsigma = zeros(Mg, 1);
	vec RinSmu_j;
	double RinSmuj_j;
	double alphag = 0;
	double sigma2e = 1;
	double a = 1;
	double sigma2mu = 1;
	double Ex2 = 0;


	while (iter <= max_iter && abs(diff) > tol)
	{
		// E step
		xxsigma = diagxtx / sigma2e + alphag*alphag*diaginSRinS + 1 / sigma2mu;
		v2 = 1. / xxsigma;

		for (int j = 0; j < Mg; j++){
			ytilde_j = ytilde - x.col(j)*mu(j);
			RinSmu_j = RinSmu - RinS.col(j)*mu(j);
			RinSmuj_j = RinSmu(j) - R(j, j)*mu(j) / hats(j);
			mu(j) = (sum(x.col(j) % (y - W*beta - ytilde_j)) / sigma2e + alphag * hatmus2(j) - alphag * alphag * RinSmuj_j / hats(j)) / xxsigma(j);
			RinSmu = RinSmu_j + RinS.col(j)*mu(j);
			ytilde = ytilde_j + x.col(j)*mu(j);
		}

		// M step
		muinSRinSmu = (mu / hats).t()*RinSmu;
		Ex2 = as_scalar(mu.t()*mu + sum(v2));

		// update sigma2e2
		sigma2e = (pow(norm(y - W*beta - ytilde, 2), 2) + sum(v2%diagxtx)) / n;

		// update sigma2mu
		sigma2mu = Ex2 / Mg;

		// update alphag
		if (fix_alphag == 0){
			alphag = as_scalar(sum(hatmus2%mu) / (muinSRinSmu + sum(v2%diaginSRinS)));
		}
		else{
			alphag = 0;
		}

		// update beta2
		beta = solve(W, y - ytilde);

		// update px
		if (px_ == 1){
			gam = as_scalar(sum(ytilde%(y - W*beta)) / (ytilde.t()*ytilde + sum(v2%diagxtx)));
			sigma2mu = gam * gam * sigma2mu;
			alphag = alphag / gam;
			mu = mu * gam;
			ytilde = ytilde * gam;
			RinSmu = RinSmu * gam;
			v2 = v2*gam*gam;
			muinSRinSmu = muinSRinSmu * gam *gam;
			Ex2 = Ex2 * gam * gam;
		}

		// lower bound
		term1_1 = as_scalar(-0.5*n*log(sigma2e) - 0.5 / sigma2e*pow(norm(y - W*beta - ytilde), 2));
		term1 = term1_1 - 0.5 / sigma2e*sum(v2%diagxtx);
		term2 = as_scalar(alphag*hatmus2.t()*mu - 0.5*alphag*alphag*muinSRinSmu - 0.5*alphag*alphag*sum(v2%diaginSRinS));
		term3 = -0.5*Mg*log(sigma2mu) - 0.5*Ex2 / sigma2mu;
		term4 = 0.5*sum(log(v2));
		lowerbound(iter - 1) = term1 + term2 + term3 + term4;

		// display the value of the lower bound if needed
		F = lowerbound(iter - 1);
		diff = F - Fold;
		if (dispF == 1)
		{
			if (fmod(iter, dispEvery) == 0)
			{
				Pr(0, 0) = iter;
				Pr(0, 1) = F;
				Pr(0, 2) = Fold;
				Pr(0, 3) = diff;
				Pr.print("***Iteration*******Fnew********Fold**********Diff***");
			}
		}
		Fold = F;

		if (diff<0)
		{
			throw std::range_error("The lower bound fails to increase");
		}

		iter = iter + 1;
	}

	// loglikelihood
	double loglikelihood = LLD(muinSRinSmu, xtx, hatmus2, inSRinS, mu, alphag, sigma2mu, Mg, sigma2e, term1_1, v2) - as_scalar(0.5*hatmu.t()*invSRS*hatmu);

	ObjSCLMM obj;

	obj.vardist_mu = mu;
	obj.alphag = alphag;
	obj.a = a;
	obj.sigma2mu = sigma2mu;
	obj.sigma2beta = sigma2beta;
	obj.sigma2y = sigma2y;
	obj.LRLB = loglikelihood;
	obj.Lq = lowerbound(span(0, iter - 2));

	return obj;
}
