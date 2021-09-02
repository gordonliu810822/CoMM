#include "RcppArmadillo.h"
#include "math.h"
#include <cassert>
#include <cmath>
#include "CoMM_S4.hpp"

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



ObjCoMM_S4 rcpparma_CoMM_S4_VB(vec &hatmur, vec &hatmu2r, mat &R, mat &R2, Options_CoMM_S4* opts, bool px){

	// check number of input arguments
	bool dispF = opts->dispF;
	int dispEvery = opts->display_gap;
	int max_iter = opts->max_iter;
	double tol = opts->epsStopLogLik;
	bool fix_alphag = opts->fix_alphag;

	double diff = datum::inf;
	double Fold = -datum::inf;

	int Mg = hatmur.n_rows;

	colvec hatmu(hatmur.begin(), hatmur.size(), false);
	colvec hatmu2(hatmu2r.begin(), hatmu2r.size(), false);

	int iter = 1;
	colvec mu = zeros(Mg, 1); // initial variational mean
	rowvec lowerbound = zeros(1, max_iter + 1);

	double term1, term2, term3, term4;
	double F = 0;
	rowvec Pr = zeros(1, 4);

	mat SIGMA = zeros(Mg, Mg);
	vec v2 = zeros(Mg, 1);

	double alphag = 0;
	double sigma2mu = 1;

	// predefined constant
	vec diagR = diagvec(R);
	vec diagR2 = diagvec(R2);
	vec Rmu = R*mu;
	vec R2mu = R2*mu;

	mat muRmu;
	mat muR2mu;
	vec xxsigma = zeros(Mg, 1);
	vec Rmu_j, R2mu_j;
	double Rmuj_j, R2muj_j;
	double LRLB;
	double a = 1;

	while (iter <= max_iter && abs(diff) > tol)
	{
		// E step
		xxsigma = a*a*diagR + alphag*alphag*diagR2 + 1 / sigma2mu;
		v2 = 1. / xxsigma;

		for (int j = 0; j < Mg; j++){
			Rmu_j = Rmu - R.col(j)*mu(j);
			R2mu_j = R2mu - R2.col(j)*mu(j);
			Rmuj_j = Rmu(j) - R(j, j)*mu(j);
			R2muj_j = R2mu(j) - R2(j, j)*mu(j);
			mu(j) = (a*hatmu(j) - a*a*Rmuj_j + alphag*hatmu2(j) - alphag * alphag * R2muj_j) / xxsigma(j);
			Rmu = Rmu_j + R.col(j)*mu(j);
			R2mu = R2mu_j + R2.col(j)*mu(j);
		}

		muRmu = mu.t()*R*mu;
		muR2mu = mu.t()*R2*mu;

		// M step
		// update sigma2mu
		sigma2mu = (sum(mu%mu) + sum(v2)) / Mg;

		// update alphag
		if (fix_alphag == 0){
			alphag = as_scalar((mu.t()*hatmu2) / (muR2mu + sum(v2%diagR2)));
		}
		else{
			alphag = 0;
		}
		if (px == 1){
			a = as_scalar((mu.t()*hatmu) / (muRmu + sum(v2%diagR)));
		}


		// lower bound
		term1 = as_scalar(a*mu.t()*hatmu - 0.5*a*a*muRmu - 0.5*a*a*sum(v2%diagR));
		term2 = as_scalar(alphag*mu.t()*hatmu2 - 0.5*alphag*alphag*muR2mu - 0.5*alphag*alphag*sum(v2%diagR2));
		term3 = -0.5*Mg*sum(log(sigma2mu)) - 0.5*sum((mu%mu + v2) / sigma2mu);
		term4 = 0.5*sum(log(v2));
		//cout << term1 << endl;
		//cout << term2 << endl;
		//cout << term3 << endl;
		//cout << term4 << endl;

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

	// linear response
	mu = mu*a;
	sigma2mu = sigma2mu*a*a;
	alphag = alphag / a;
	SIGMA = inv_sympd(R + alphag*alphag*R2 + 1 / sigma2mu*eye(Mg, Mg));
	term1 = as_scalar(mu.t()*hatmu - 0.5*mu.t()*R*mu - 0.5*trace(SIGMA*R));
	term2 = as_scalar(alphag*mu.t()*hatmu2 - 0.5*alphag*alphag*mu.t()*R2*mu - 0.5*alphag*alphag*trace(SIGMA*R2));
	term3 = -0.5*Mg*sum(log(sigma2mu)) - 0.5*sum((mu%mu + diagvec(SIGMA)) / sigma2mu);
	term4 = sum(log(diagvec(chol(SIGMA))));
	LRLB = term1 + term2 + term3 + term4;

	ObjCoMM_S4 obj;

	obj.vardist_mu = mu;
	obj.alphag = alphag;
	obj.sigma2mu = sigma2mu;
	obj.LRLB = LRLB;
	obj.Lq = lowerbound(span(0, iter - 2));

	return obj;
}