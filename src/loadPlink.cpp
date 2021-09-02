#include <stdio.h>
#include <bitset>
#include <math.h>
#include <iostream>
#include "plinkfun.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;

void loadplink(std::string stringname) {
  string famfile = stringname;
  famfile += ".fam";
  string bimfile = stringname;
  bimfile += ".bim";
  int N =  getLineNum(famfile);
  int P =  getLineNum(bimfile);
  clock_t t1 = clock();
  //Chroms chroms(bimfile, P);
  //int phenotype_pos = 6;
  //vec y = read_phenotypes(famfile, N, phenotype_pos);
//  unsigned* X = new unsigned[ N * P];
  arma::Mat<unsigned> Xdata(N,P);
  t1 = clock();
  readPlink(stringname,N, P, Xdata.memptr());
 // arma::Mat<unsigned>* Xdata = new arma::Mat<unsigned>(X, N, P, false,false);
  cout<<"Finish Reading Plink file, time elapsed:"<< (clock() - t1)*1.0/CLOCKS_PER_SEC << endl;
  cout <<"Sample Size = " <<  N << " SNP Number = " << P << endl;
 // delete[] X;
}
