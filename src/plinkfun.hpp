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

void readPlink(string stringname, long int N, long int P, unsigned* X);
void readPlink2(string stringname,int N, int P, char* X);
int getLineNum(string filename);
arma::mat getSubMat(char* X, int N , int P, uvec indics, double* sub_matrix_double);
void ReadPlinkBimFile(string stringname, IntegerVector A1, IntegerVector A2, CharacterVector rsname,
                      IntegerVector chr, IntegerVector bp, NumericVector morgan, int P);
void ReadPlinkFamFile(std::string stringname, CharacterVector FID, CharacterVector IID, IntegerVector sex,
                      NumericVector pheno, int N);
void ReadPlinkFamFile2(std::string stringname, CharacterVector FID, CharacterVector IID,
                       NumericVector pheno, int nrows, int whCol);
//void getFourGentype(unsigned* geno, std::bitset<8> bits);
CharacterVector charv_subset(CharacterVector x, uvec idx);

#endif /* plinkfun_hpp */
