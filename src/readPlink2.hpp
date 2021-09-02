//  readPlink.cpp
//  ReadPlinkGit

#ifndef readPlink2_hpp
#define readPlink2_hpp

#include <RcppArmadillo.h>
#include<iostream>
#include <stdio.h>
#include "plinkfun.hpp"


using namespace std;
using namespace arma;

struct ObjXY{
    arma::Mat<unsigned> X;
	arma::colvec Y;
};

ObjXY ReadDataFromFile(string stringname);


#endif
