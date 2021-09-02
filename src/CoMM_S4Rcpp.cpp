#include "RcppArmadillo.h"
#include "CoMM_S4.hpp"
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
List CoMM_S4(arma::vec &hatmur, arma::vec &hatmu2r, arma::mat &Rr, arma::mat &R2r, SEXP opts = R_NilValue, bool px = 1) {

	int Mg = hatmur.n_rows;

	colvec hatmu(hatmur.begin(), hatmur.size(), false);
	colvec hatmu2(hatmu2r.begin(), hatmu2r.size(), false);
	mat R(Rr.begin(), Mg, Mg, false);
	mat R2(R2r.begin(), Mg, Mg, false);

	Options_CoMM_S4* lp_opt = NULL;
	if (!Rf_isNull(opts)){
		Rcpp::List opt(opts);
		lp_opt = new Options_CoMM_S4(opt["max_iter"], opt["dispF"], opt["display_gap"], opt["epsStopLogLik"], opt["fix_alphag"]);
	}

	if (Rf_isNull(opts)){
		lp_opt = new Options_CoMM_S4();
	}

	ObjCoMM_S4 obj = rcpparma_CoMM_S4_VB(hatmu, hatmu2, R, R2, lp_opt, px);


	List ret;
	ret["vardist_mu"] = Rcpp::wrap(obj.vardist_mu);
	ret["sigma2mu"] = Rcpp::wrap(obj.sigma2mu);
	ret["alphag"] = Rcpp::wrap(obj.alphag);
	ret["Lq"] = Rcpp::wrap(obj.Lq);
	ret["LRLB"] = Rcpp::wrap(obj.LRLB);
	return ret;
}