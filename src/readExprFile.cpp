#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include "plinkfun.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace Rcpp;
using namespace arma;

Rcpp::List getExprFile2(std::string filename, char delimiter, int ncols, int nrows, int ncols_omit){
  std::ifstream myfile (filename.c_str());
  mat expr(nrows -1, ncols-ncols_omit);
  std::string line;
  CharacterVector genetype1(nrows - 1), genetype2(nrows - 1), targetID(nrows - 1);
  vec lower(nrows - 1), upper(nrows - 1), chr_expr(nrows - 1);

  if (myfile.is_open()){
    getline(myfile, line);
  }
  clock_t t1 = clock();

  int nrow_ind = 0;
  vector <string> tmp;

  while (nrow_ind < nrows - 1){
    if (nrow_ind % 5000 == 0 && nrow_ind != 0){
      cout << nrow_ind << "-th Gene" << ",";
      cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
    }

    getline(myfile, line);

    boost::split( tmp, line, boost::is_any_of(" \t *"));

    genetype1(nrow_ind) = tmp[2];
    genetype2(nrow_ind) = tmp[3];
    targetID(nrow_ind) = tmp[4];
    lower(nrow_ind) = atof(tmp[0].c_str());
    upper(nrow_ind) = atof(tmp[1].c_str());
    chr_expr(nrow_ind) = atof(tmp[5].c_str());


    for (int j = 0; j < ncols-ncols_omit; j ++){
      if (tmp[j+6].compare("NA") == 0){
        expr(nrow_ind,j) = log(-10);
      }
      else {
        expr(nrow_ind,j) = atof(tmp[j+6].c_str());
      }
    }

    nrow_ind ++ ;

  }

  List output = List::create(Rcpp::Named("genetype1") = genetype1,
                             Rcpp::Named("genetype2") = genetype2,
                             Rcpp::Named("targetID") = targetID,
                             Rcpp::Named("lower") = lower,
                             Rcpp::Named("upper") = upper,
                             Rcpp::Named("chr_expr") = chr_expr,
                             Rcpp::Named("expr") = expr);
  return output;
}



