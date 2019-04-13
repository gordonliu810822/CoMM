#include "RcppArmadillo.h"
#include "SCLMM.hpp"
#include <Rcpp.h>

using namespace Rcpp;
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
//



// [[Rcpp::export]]
List CoMM_S2(arma::mat xr, arma::vec yr, arma::vec Wr, arma::vec &hatmur, arma::vec &hatsr, arma::mat &Rr, SEXP opts = R_NilValue, bool px = 1) {

	int Mg = hatmur.n_rows;
	int n = xr.n_rows;
	int q = Wr.n_cols;

	mat x(xr.begin(), n, Mg, false);   // pointer
	mat W(Wr.begin(), n, q, false);   // pointer
	vec y(yr.begin(), yr.size(), false);
	vec hatmu(hatmur.begin(), hatmur.size(), false);
	vec hats(hatsr.begin(), hatsr.size(), false);
	mat R(Rr.begin(), Mg, Mg, false);

	Options_SCLMM* lp_opt = NULL;
	if (!Rf_isNull(opts)){
		Rcpp::List opt(opts);
		lp_opt = new Options_SCLMM(opt["max_iter"], opt["dispF"], opt["display_gap"], opt["epsStopLogLik"], opt["fix_alphag"]);
	}

	if (Rf_isNull(opts)){
		lp_opt = new Options_SCLMM();
	}

	bool px_ = px;
	ObjSCLMM obj = rcpparma_SCLMM_IS(x, y, W, hatmu, hats, R, lp_opt, px_);

	List ret;
	ret["vardist_mu"] = Rcpp::wrap(obj.vardist_mu);
	ret["sigma2mu"] = Rcpp::wrap(obj.sigma2mu);
	ret["alphag"] = Rcpp::wrap(obj.alphag);
	ret["sigma2beta"] = Rcpp::wrap(obj.sigma2beta);
	ret["sigma2y"] = Rcpp::wrap(obj.sigma2y);
	ret["LRLB"] = Rcpp::wrap(obj.LRLB);
	ret["Lq"] = Rcpp::wrap(obj.Lq);
	return ret;

}

