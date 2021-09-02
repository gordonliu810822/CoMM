#ifndef CoMM_S4_hpp
#define CoMM_S4_hpp

#include <RcppArmadillo.h>
#include<iostream>
#include<fstream>
#include <stdio.h>

using namespace std;
using namespace arma;

class Options_CoMM_S4{
public:
	// Constructor definition
	// The complier deciedes which constructor to be called depending on 
	// the number of argument present with the object
	Options_CoMM_S4(){
		this->max_iter = 1e5;
		this->dispF = 1;
		this->display_gap = 10;
		this->epsStopLogLik = 1e-5;
		this->fix_alphag = 0;
	}
	Options_CoMM_S4(int max_iter, bool dispF, int display_gap, double epsStopLogLik, bool fix_alphag){
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


struct ObjCoMM_S4{
	colvec vardist_mu;
	double alphag;
	colvec sigma2mu;
	double sigma2e2;
	colvec beta;
	double LRLB;
	rowvec Lq;
};



ObjCoMM_S4 rcpparma_CoMM_S4_VB(vec &hatmur, vec &hatmur2, mat &R, mat &R2, Options_CoMM_S4* opts, bool px);
#endif
