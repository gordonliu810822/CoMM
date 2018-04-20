//
//  plinkfun.hpp
//  PlinkRead
//
//  Created by DaiMingwei on 16/10/6.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#ifndef plinkfun_hpp
#define plinkfun_hpp

#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <bitset>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace Rcpp;
using namespace arma;


void getFourGentype(unsigned* geno, std::bitset<8> bits);
void readPlink(string stringname, long int N, long int P, unsigned* X);
int getLineNum(string filename);
void ReadPlinkBimFile(string stringname, IntegerVector A1, IntegerVector A2, CharacterVector rsname,
	IntegerVector chr, IntegerVector bp, NumericVector morgan, int P);
void ReadPlinkFamFile(std::string stringname, CharacterVector FID, CharacterVector IID, IntegerVector sex,
	NumericVector pheno, int N);
//CharacterVector charv_subset(CharacterVector x, uvec idx);
List match_SNPs(string stringname1, string stringname2);
void ReadPlinkFamFile2(std::string stringname, CharacterVector FID, CharacterVector IID,
	NumericVector pheno, int nrows, int whCol);


#endif /* plinkfun_hpp */
