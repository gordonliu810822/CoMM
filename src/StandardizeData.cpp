#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;


arma::mat StandardX(arma::Mat<unsigned>* X){
	uword size1 = X->n_rows;
	uword size2 = X->n_cols;

	arma::mat standardX;
	standardX.set_size(size1, size2);
	X->replace(3, 0); //replace the missing value indicated by 3
	vec meanX(size2);
	vec sqrtsum(size2);
	for (unsigned i = 0; i < size2; i++) { //calculate the mean of the vector and sqrt sum
		meanX[i] = sum(X->col(i))*1.0 / size1;
		vec v_i = conv_to<vec>::from(X->col(i));
		v_i = v_i - meanX[i];
		mat pd = v_i.t() * v_i;
		sqrtsum[i] = sqrt(pd.at(0) / size1);
	}

	for (unsigned k = 0; k < size2; k++) {//standardization
		for (unsigned j = 0; j < size1; j++){
			standardX(j, k) = ((*X)(j, k) - meanX[k]) / sqrtsum[k];
		}

	}
	return standardX;
}