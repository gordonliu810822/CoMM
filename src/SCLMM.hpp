#ifndef SCLMM_hpp
#define SCLMM_hpp

#include <RcppArmadillo.h>
#include<iostream>
#include<fstream>
#include <stdio.h>

using namespace std;
using namespace arma;

class Options_SCLMM{
public:
	// Constructor definition
	// The complier deciedes which constructor to be called depending on 
	// the number of argument present with the object
	Options_SCLMM(){
		this->max_iter = 1e5;
		this->dispF = 1;
		this->display_gap = 10;
		this->epsStopLogLik = 1e-5;
		this->fix_alphag = 0;
	}
	Options_SCLMM(int max_iter, bool dispF, int display_gap, double epsStopLogLik, bool fix_alphag){
		this->max_iter = max_iter;
		this->dispF = dispF;
		this->display_gap = display_gap;
		this->epsStopLogLik = epsStopLogLik;
		this->fix_alphag = fix_alphag;
	}

	int max_iter;
	bool dispF;
	int display_gap;
	double epsStopLogLik;
	bool fix_alphag;
};


struct ObjSCLMM{
	colvec vardist_mu;
	double alphag;
	double a;
	double sigma2mu;
	double sigma2e2;
	double sigma2beta;
	double sigma2y;
	colvec beta;
	double LRLB;
	rowvec Lq;
};


ObjSCLMM rcpparma_SCLMM_IS(mat &xr, vec &yr, mat &Wr, vec &hatmur, vec &hatsr, mat &R, Options_SCLMM* opts, bool px);

#endif
