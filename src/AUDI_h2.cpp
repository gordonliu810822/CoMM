#include <RcppArmadillo.h>
//#include <Rcpp.h>
//#include <omp.h>
#include <stdio.h>
#include <bitset>
#include <math.h>

#include "data_loader.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;

template <class T> const T& max3(const T& a, const T& b) {
	return (a < b) ? b : a;     // or: return comp(a,b)?b:a; for version (2)
}

template <class T> const T& min3(const T& a, const T& b) {
	return (a > b) ? b : a;     // or: return comp(a,b)?b:a; for version (2)
}

// [[Rcpp::export]]
Rcpp::List mmhe(arma::vec y, arma::mat X, arma::mat K){
	uword n_subj = y.n_elem;
	uword n_cov = X.n_cols;
	double trK = trace(K), trKK = sum( sum(K % K));
	rowvec yK = trans(y)*K;
	mat XK = X.t()*K;

	mat XX = X.t()*X, Z = X*XX.i();
	rowvec yZ = trans(y)*Z;
	double yPy = sum(y % y) - sum(yZ*X.t()*y);

	rowvec yZXK = yZ*XK;
	mat XKZ = XK*Z;

	double yPKPy = sum(yK*y)- 2 * sum(yZXK*y) + sum(yZXK*X*yZ.t());
	double trPK = trK - trace(XKZ);
	double trPKPK = trKK-2*trace(XX.i()*XK*XK.t())+trace(XKZ*XKZ);

	mat S = zeros<mat>(2, 2);
	S(0,0) = trPKPK;
	S(0,1) = trPK;
	S(1,0) = trPK;
	S(1,1) = n_subj-n_cov;
	vec q = zeros<vec>(2);
	q(0) = yPKPy;
	q(1) = yPy;

	vec Vc = S.i()*q;

	Vc.elem( find(Vc < 0) ).zeros();
	double s = trPKPK-trPK*trPK/((double)n_subj-(double)n_cov);
	double h2 = max3(min3(Vc(0)/sum(Vc),1.0),0.0);
	double se = sqrt(2/s);

	List out = List::create(Rcpp::Named("Vc") = Vc,
		Rcpp::Named("s") = s,
		Rcpp::Named("h2") = h2,
		Rcpp::Named("se") = se,
		Rcpp::Named("diff") = sum(S,1) - q);


	return out;
}



// [[Rcpp::export]]
Rcpp::List AUDI_h2_run(std::string stringname1, std::string stringname2, std::string stringname3, std::string stringname4,
	std::string stringname5, int whCol, int bw){ //int normalize_option = 1, int pred_option = 0){//, char* A21, char* A22){
	// normalize_option: 1. normalize each separately, 2. normalize both plink files together
	// match SNPs in file 1 and file 2 GWAS (common SNPs in x1 and x2 in columns)
	// plink file 1: stringname1; plink file 2: stringname2; expression file: stringname3
	// covariates file for file 1: stringname4; covariates file for file 2: stringname5
	// pred_option :0 (no calculation for prediction) 1 (calcuation for prediction)

	List tmp = dataLoader(stringname1, stringname2, stringname3, stringname4, stringname5, whCol);
	Mat<unsigned> X1 = tmp["X1"], X2 = tmp["X2"];
	vec z = tmp["y"];
	mat expr = tmp["expr_used"];
	mat w1 = tmp["covar1"], w2 = tmp["covar2"];
	uvec chr_4use_r = tmp["chr_4use_r"];
	uvec bp_4use_r = tmp["bp_4use_r"];
	vec lower = tmp["lower"], upper = tmp["upper"], chr_expr = tmp["chr_expr"];
	uword Ngene = lower.size();

	double trXmuS = 0, trXXmuS = 0;
	uvec idx;

	// Choose SNPs to calculate K;
	mat K1 = zeros<mat>(z.n_elem, z.n_elem);
	mat K2 = zeros<mat>(z.n_elem, z.n_elem); //for all snps h2

	cout << "## Start calculate K ... "  << endl;
	clock_t t0 = clock();

	//mat X2_all = zeros<mat>(z.n_elem, 0);
	// mat X2_all2 = zeros<mat>(z.n_elem, 0);
	// mat X1_all = zeros<mat>(X1.n_rows, 0);
	uvec idx_all = zeros<uvec>(0);
	//vec idx_start = zeros<vec>(Ngene), idx_end = zeros<vec>(Ngene);
	uvec idx_gene = zeros<uvec>(Ngene);
	int pX2 = 0;
	//cout << "n2: " << X1.n_rows << endl;

	for (uword g = 0; g < Ngene; g++){

		idx = find(conv_to<vec>::from(bp_4use_r) < upper(g) + bw && conv_to<vec>::from(bp_4use_r) > lower(g) - bw
			&& conv_to<vec>::from(chr_4use_r) == chr_expr(g));

		if (g % 1000 == 0 && g != 0){
			cout << g + 1 << "-th Gene starts working ...";
			cout << "Elapsed time is " << (clock() - t0)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
		}

		if (idx.is_empty() == false){
			if (idx.n_elem > 1){

				mat X1tmp = conv_to<mat>::from((&X1)->cols(idx));
				mat X2tmp = conv_to<mat>::from((&X2)->cols(idx));

				uvec idx1(1); idx1(0) = g;
				vec y = trans(expr.rows(idx1));
				y.replace(datum::nan, 9999);

				uvec idx2 = find(y != 9999);
				X1tmp = ((&X1tmp)->rows(idx2));
				y = y(idx2);

				mat w1tmp = ((&w1)->rows(idx2));

				rowvec meanX1tmp = mean(X1tmp, 0);
				rowvec sdX1tmp = stddev(X1tmp, 0, 0); // see manual
				X1tmp = (X1tmp - repmat(meanX1tmp, X1tmp.n_rows, 1)) / repmat(sdX1tmp, X1tmp.n_rows, 1) / sqrt(X1tmp.n_cols);

				//cout << "break 1..." << endl;

				/*List fm0 = lmm_pxem(y, w1tmp, X1tmp, 100);

				double sigma2y = fm0["sigma2y"], sigma2beta = fm0["sigma2beta"];
				vec beta0 = fm0["beta0"], loglik = fm0["loglik"];
				mat Sigb = fm0["Sigb"];
				vec mub = fm0["mub"];
				double h2tmp = sigma2beta / (sigma2beta + sigma2y);*/

				List fm0 = mmhe(y, w1tmp, X1tmp*X1tmp.t());
				double h2tmp = fm0["h2"];

				//mat X2tmp = conv_to<mat>::from((&X2) ->cols(idx));
				if (h2tmp > 0.01){

					idx_gene(g) = 1;
					idx_all = join_cols(idx_all, idx);
					/*if (g != 0){
						idx_start(g) = idx_end(g - 1) + idx_gene(g-1);
						// if g-1 gene used, start at idx_end(g - 1) + 1, else start at idx_end(g - 1).
						idx_end(g) = idx_start(g) + idx.n_elem - 1;
					}
					else{
						idx_start(g) = 0;
						idx_end(g) = idx.n_elem - 1;
					}	*/

					rowvec meanX2tmp = mean(X2tmp, 0);
					rowvec sdX2tmp = stddev(X2tmp, 0, 0); // see manual
					X2tmp = (X2tmp - repmat(meanX2tmp, X2tmp.n_rows, 1)) / repmat(sdX2tmp, X2tmp.n_rows, 1);

					//X2_all = join_rows(X2_all, X2tmp / sqrt(X2tmp.n_cols));
					// X2_all2 = join_rows(X2_all2, X2tmp);

					mat tmp = X2tmp*X2tmp.t();
					K2 = tmp + K2;
					K1 = tmp/X2tmp.n_cols + K1;
					
					pX2 = pX2 + (int)X2tmp.n_cols;
					//X1_all = join_rows(X1_all, X1tmp);
				}
				/*else{
					if (g != 0){
						idx_start(g) = idx_start(g - 1);
						idx_end(g) = idx_end(g - 1);
					}
					else{
						idx_start(g) = 0;
						idx_end(g) = 0;
					}

				}*/
			}
		}
	}
	cout << "sum of idx_gene: " << sum(idx_gene) << endl;
	/*List out = List::create(Rcpp::Named("X2_all") = X2_all,
		Rcpp::Named("X2_all2") = X2_all2,
		Rcpp::Named("idx_all") = idx_all,
		Rcpp::Named("idx_gene") = idx_gene);*/

	// int nAllCols = X2_all.n_cols;
    // uvec idxall = linspace<uvec>(1, nAllCols - 1, nAllCols - 1);
	// X2_all = X2_all.cols(idxall);

	// List fm0 = lmm_pxem(y, w1tmp, X1tmp, 100);

	//X2_all = X2_all / sqrt(nAllCols - 1);
	//cout << "Dim of X2_all: " << X2_all.n_rows << "*" << X2_all.n_cols << endl;
	cout << "The dimension of X2_all is " << X2.n_rows << "*" << pX2 << "," << idx_all.n_elem << endl;
    //K = X2_all *X2_all.t();

	List fmall = mmhe(z, w2, K2/pX2);
	double h20 = fmall["h2"];
	vec Vc = fmall["Vc"];
	//List fmall = lmm_pxem(z, w2, X2_all2, 200);
	//double sigma2z0 = fmall["sigma2y"], sigma2alpha = fmall["sigma2beta"];
	//double h20 = sigma2alpha / (sigma2alpha + sigma2z0);
	/*mat X2tmp = conv_to<mat>::from(X2);
	rowvec meanX2tmp = mean(X2tmp, 0);
	rowvec sdX2tmp = stddev(X2tmp, 0, 0); // see manual
	X2tmp = (X2 - repmat(meanX2tmp, X2tmp.n_rows, 1))/ repmat(sdX2tmp, X2tmp.n_rows, 1) / sqrt(X2tmp.n_cols);
	K = X2tmp * X2tmp.t();*/



	double trK = trace(K1);
	cout << " Elapsed time is " << (clock() - t0)*1.0 / CLOCKS_PER_SEC << " sec" << endl;


	for (uword g = 0; g < Ngene; g++){
		idx = find(conv_to<vec>::from(bp_4use_r) < upper(g) + bw && conv_to<vec>::from(bp_4use_r) > lower(g) - bw
			&& conv_to<vec>::from(chr_4use_r) == chr_expr(g));

		if ( g % 1000 == 0 && g != 0){
					cout << g + 1 << "-th Gene starts working ..." ;
					cout << "Elapsed time is " << (clock() - t0)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
		}

		if (idx_gene(g) == 1){
			mat X1tmp = conv_to<mat>::from((&X1)->cols(idx));
			mat X2tmp = conv_to<mat>::from((&X2)->cols(idx));

			uvec idx1(1); idx1(0) = g;
			vec y = trans(expr.rows(idx1));
			y.replace(datum::nan, 9999);

			uvec idx2 = find(y != 9999);
			X1tmp = ((&X1tmp)->rows(idx2));

			mat w1tmp = ((&w1)->rows(idx2));

			y = y(idx2);

			rowvec meanX1tmp = mean(X1tmp, 0);
			rowvec sdX1tmp = stddev(X1tmp, 0, 0); // see manual
			X1tmp = (X1tmp - repmat(meanX1tmp, X1tmp.n_rows, 1)) / repmat(sdX1tmp, X1tmp.n_rows, 1) / sqrt(X1tmp.n_cols);

			rowvec meanX2tmp = mean(X2tmp, 0);
			rowvec sdX2tmp = stddev(X2tmp, 0, 0); // see manual
			X2tmp = (X2tmp - repmat(meanX2tmp, X2tmp.n_rows, 1)) / repmat(sdX2tmp, X2tmp.n_rows, 1) / sqrt(X2tmp.n_cols);

			List fm0 = lmm_pxem(y, w1tmp, X1tmp, 100);

			//double sigma2y = fm0["sigma2y"], sigma2beta = fm0["sigma2beta"];
			//vec beta0 = fm0["beta0"], loglik = fm0["loglik"];
			mat Sigb = fm0["Sigb"];
			vec mub = fm0["mub"];

			// Calculate MoM of h2
			vec x2mu = X2tmp*mub;
			mat x2tx2 = X2tmp.t()*X2tmp;
			trXmuS = trXmuS + sum(x2mu % x2mu) + trace(x2tx2*Sigb);
			trXXmuS = trXXmuS + sum(x2mu.t()*K1*x2mu) + trace(X2tmp.t() *K1 * X2tmp *Sigb);

		}
	}
	//delete[] X2_all;

	// Calculate MoM of h2
	mat S = zeros<mat>(2, 2);
	S(0,0) = trXmuS;
	S(0,1) = z.n_elem;
	S(1,0) = trXXmuS;
	S(1,1) = trK;
	vec q = zeros<vec>(2);
	q(0) = sum(z % z);
	//vec x2tz = X2_all.t()*z;
	//q(1) = sum(x2tz % x2tz);
	mat kernel = z.t() *K1 *z;
	q(1) = kernel(0,0);

	vec varcomp = solve( S, q ) ;
	double sigma2a = varcomp[0], sigma2z = varcomp[1];
	double m = sum(idx_gene);
	double h22 = m /(m+ sigma2z/sigma2a);

	List out = List::create(Rcpp::Named("varcomp") = varcomp,
		Rcpp::Named("Vc") = Vc,
		//Rcpp::Named("sigma2a") = sigma2a,
		Rcpp::Named("h22") = h22,
		Rcpp::Named("h20") = h20,
		Rcpp::Named("z") = z,
		Rcpp::Named("idx_all") = idx_all,
		Rcpp::Named("K1") = K1,
		Rcpp::Named("K2") = K2,
		Rcpp::Named("S") = S,
		Rcpp::Named("q") = q,
		Rcpp::Named("idx_gene") = idx_gene);
		//Rcpp::Named("X2_all") = X2_all);

	return out;
}
