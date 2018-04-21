#include <RcppArmadillo.h>
#include <Rcpp.h>
//#include <omp.h>
#include <stdio.h>
#include <bitset>
#include <math.h>


#include "data_loader.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;

//' @title
//' CoMM
//' @description
//' CoMM to dissecting genetic contributions to complex traits by leveraging regulatory information.
//'
//' @param stringname1  prefix for eQTL genotype file with plink format (bim, bed).
//' @param stringname2  prefix for GWAS genotype and phenotype file with plink format (bim, bed, fam).
//' @param stringname3  gene expression file with full name.
//' @param stringname4  covariates file for eQTL data.
//' @param stringname5  covariates file for GWAS data, e.g. top 10 PCs.
//' @param whCol  specify which phenotype is used in fam. For example, when whCol = 2, the seven-th column of fam file will be used as phenotype.
//' @param bw  the number of downstream and upstream SNPs that are considered as cis-SNP within a gene.
//'
//' @return List of model parameters
//'
//' @examples
//' ##Working with no summary statistics, no covariates and options
//' data(CoMM)
//' file1 = "1000G.EUR.QC.1";
//' file2 = "NFBC_filter_mph10";
//' file3 = "Geuvadis_gene_expression_qn.txt";
//' file4 = "";
//' file5 = "pc5_NFBC_filter_mph10.txt";
//' whichPheno = 1;
//' bw = 500000;
//'
//' fm = AUDI_testing_run(file1,file2,file3, file4,file5, whichPheno, bw);
//' data(IGESSDB)
//' fit <- iGess(Model$X, Model$y, NULL, Model$AIPVal)
//'
//' @details
//' \code{CoMM} fits the CoMM model. It requires to provide plink binary eQTL genotype file (bim, bed)
//' the GWAS plink binary file (bim, bed, fam), gene expression file for eQTL.
//' @export
// [[Rcpp::export]]
Rcpp::List CoMM_testing_run(std::string stringname1, std::string stringname2, std::string stringname3, std::string stringname4,
	std::string stringname5, int whCol, int bw){ //int normalize_option = 1, int pred_option = 0){//, char* A21, char* A22){
	// normalize_option: 1. normalize each separately, 2. normalize both plink files together
	// match SNPs in file 1 and file 2 GWAS (common SNPs in x1 and x2 in columns)
	// plink file 1: stringname1; plink file 2: stringname2; expression file: stringname3
	// covariates file for file 1: stringname4; covariates file for file 2: stringname5
	// pred_option :0 (no calculation for prediction) 1 (calcuation for prediction)

	List tmp = dataLoader(stringname1, stringname2, stringname3, stringname4, stringname5, whCol);
	Mat<unsigned> X1 = tmp["X1"], X2 = tmp["X2"];
	vec z = tmp["y"];
	mat w1 = tmp["covar1"], w2 = tmp["covar2"];
	mat expr = tmp["expr_used"];
	CharacterVector rsname_4use_r = tmp["rsname_4use_r"];
	uvec chr_4use_r = tmp["chr_4use_r"];
	uvec bp_4use_r = tmp["bp_4use_r"];
	CharacterVector genetype1 = tmp["genetype1"], genetype2 = tmp["genetype2"], targetID = tmp["targetID"];
	vec lower = tmp["lower"], upper = tmp["upper"], chr_expr = tmp["chr_expr"];

	uword Ngene = lower.size();
	umat snp_info = zeros<umat>(chr_4use_r.n_elem, 2);
	snp_info.col(0) = chr_4use_r;
	snp_info.col(1) = bp_4use_r;

	mat expr_info = zeros<mat>(lower.n_elem, 3);
	expr_info.col(0) = chr_expr;
	expr_info.col(1) = lower;
	expr_info.col(2) = upper;

	uvec idx;

	//mat y_pred = zeros<mat>(1, z.n_elem), y_se = zeros<mat>(1, z.n_elem);

	mat out_param = zeros<mat>(Ngene, 7);
	uvec idx_all = zeros<uvec>(0);


	// conduct testing

	clock_t t1 = clock();
	for (uword g = 0; g < Ngene; g++){

		idx = find(conv_to<vec>::from(bp_4use_r) < upper(g) + bw && conv_to<vec>::from(bp_4use_r) > lower(g) - bw
			&& conv_to<vec>::from(chr_4use_r) == chr_expr(g));


		out_param(g, 5) = idx.n_elem;

		if ( g % 100 == 0 && g != 0){
					cout << g + 1 << "-th Gene starts working ..." ;
					cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
		}

		if (idx.is_empty() == false){
			if (idx.n_elem > 1){
				idx_all = join_cols(idx_all, idx);

				mat X1tmp = conv_to<mat>::from((&X1) ->cols(idx));
				mat X2tmp = conv_to<mat>::from((&X2) ->cols(idx));

				uvec idx1(1); idx1(0) = g;
				vec y = trans(expr.rows(idx1));
				y.replace(datum::nan, 9999);

				uvec idx2 = find(y != 9999);
				X1tmp = ((&X1tmp) -> rows(idx2));


				mat w1tmp = ((&w1) -> rows(idx2));

				y = y(idx2);

				//if (normalize_option == 1){
				rowvec meanX1tmp = mean(X1tmp, 0);
				rowvec sdX1tmp = stddev(X1tmp, 0, 0); // see manual
				X1tmp = (X1tmp - repmat(meanX1tmp, X1tmp.n_rows, 1))/ repmat(sdX1tmp, X1tmp.n_rows, 1) / sqrt(X1tmp.n_cols);

				//cout << "Dim of X1tmp: " << X1tmp.n_rows <<"*" <<X1tmp.n_cols << endl;

				rowvec meanX2tmp = mean(X2tmp, 0);
				rowvec sdX2tmp = stddev(X2tmp, 0, 0); // see manual
				X2tmp = (X2tmp - repmat(meanX2tmp, X2tmp.n_rows, 1))/ repmat(sdX2tmp, X2tmp.n_rows, 1) / sqrt   (X2tmp.n_cols);

				rowvec X1row = X1tmp.row(0);
				X1tmp = X1tmp.cols(find_finite(X1row));
				//cout << g + 1 << "-th Gene: dim of X1 tmp " << X1tmp.n_rows << "-by-" << X1tmp.n_cols;
				//X1out = X1tmp;

				X2tmp = X2tmp.cols(find_finite(X1row));
				//cout << "Dim of X2tmp: " << X2tmp.n_rows <<"*" <<X2tmp.n_cols << endl;
				//}
				/*else if ( normalize_option == 2 ){
					mat X12 = join_cols(X1tmp, X2tmp);

					rowvec meanXtmp = mean(X12, 0);
					rowvec sdXtmp = stddev(X12, 0, 0); // see manual
					X12 = (X12 - repmat(meanXtmp, X12.n_rows, 1))/ repmat(sdXtmp, X12.n_rows, 1) / sqrt(X12.n_cols);

					uvec idxtmp1 = linspace<uvec>(0, X1tmp.n_rows - 1, X1tmp.n_rows);
					uvec idxtmp2 = linspace<uvec>(X1tmp.n_rows, X1tmp.n_rows + X2tmp.n_rows - 1, X2tmp.n_rows);

					X1tmp = X12.rows(idxtmp1);
					X2tmp = X12.rows(idxtmp2);
				}*/

				// initialize by linear mixed model for sigma2y and sigma2beta
				// use this initialization to calculate MoM estimates of h2
				List fm0 = lmm_pxem(y, w1tmp, X1tmp, 100);

				//cout <<"break 1 ..." <<endl;
				double sigma2y = fm0["sigma2y"], sigma2beta = fm0["sigma2beta"];
				vec beta0 = fm0["beta0"], loglik = fm0["loglik"];
				mat Sigb = fm0["Sigb"];
				vec mub = fm0["mub"];

				out_param(g, 6) = 1/(1+sigma2y/sigma2beta); // gene-wise h2
				/*if (pred_option == 1){
					// calculate prediction mean and variance
					y_pred = join_cols(y_pred, trans(X2tmp*mub));
					mat kernel = X2tmp * Sigb * X2tmp.t();
					vec se2 = kernel.diag() + sigma2y;

					y_se = join_cols(y_se, trans(sqrt(se2)));
				}*/

				// fit model under Ha and H0
				List fmHa = CoMM_covar_pxem(y, z, X1tmp, X2tmp, w1tmp, w2, sigma2beta, sigma2y, beta0, 0, 1e-5, 1000);
				List fmHo = CoMM_covar_pxem(y, z, X1tmp, X2tmp, w1tmp, w2, sigma2beta, sigma2y, beta0, 1, 1e-5, 1000);
				//cout <<"break 2 ..." <<endl;

				double alphatmp = fmHa["alpha"];
				vec loglikHa = fmHa["loglik"], loglikHo = fmHo["loglik"];
				sigma2y = fmHa["sigma2y"]; sigma2beta = fmHa["sigma2beta"];
				double sigma2z = fmHa["sigma2z"];

				out_param(g, 0) = sigma2beta;
				out_param(g, 1) = sigma2y;
				out_param(g, 2) = sigma2z;
				out_param(g, 3) = 2 * (max(loglikHa) - max(loglikHo));
				out_param(g, 4) = alphatmp;

			}
			else{
				out_param(g, 3) = -99;
				out_param(g, 4) = -99;
			}
		}

		else{
			out_param(g, 3) = -99;
			out_param(g, 4) = -99;
		}

	}

	/*// Calculate MoM of h2
	mat S = zeros<mat>(2, 2);
	S(0,0) = trXmuS;
	S(0,1) = z.n_elem;
	S(1,0) = trXXmuS;
	S(1,1) = trK;
	vec q = zeros<vec>(2);
	q(0) = sum(z % z);
	mat kernel = z.t() *K *z;
	q(1) = kernel(0,0);

	vec varcomp = solve( S, q ) ;
	double sigma2a = varcomp[0], sigma2z = varcomp[1];
	double h22 = 1/(1+ sigma2z/sigma2a);*/

	/*if (pred_option == 1){
		int nActiveGene = y_pred.n_rows;
		uvec idxpred = linspace<uvec>(1, nActiveGene - 1, nActiveGene - 1);
		y_pred = y_pred.rows(idxpred);
		y_se = y_se.rows(idxpred);
	}*/

	List out = List::create(//Rcpp::Named("Ngene") = Ngene,
		//Rcpp::Named("chr_expr") = chr_expr,
		//Rcpp::Named("sigma2yout") = sigma2yout,
		//Rcpp::Named("sigma2betaout") = sigma2betaout,
		//Rcpp::Named("beta0out") = beta0out,
		//Rcpp::Named("loglikout") = loglikout,idx_all
		Rcpp::Named("idx_all") = idx_all,
		Rcpp::Named("X1") = X1,
		Rcpp::Named("X2") = X2,
		Rcpp::Named("z") = z,
		Rcpp::Named("covar1") = w1,
		Rcpp::Named("covar2") = w2,
		Rcpp::Named("expr") = expr,
		Rcpp::Named("rsname_4use_r") = rsname_4use_r,
		//Rcpp::Named("chr_4use_r") = chr_4use_r,
		//Rcpp::Named("bp_4use_r") = bp_4use_r,
		Rcpp::Named("snp_info_4use_r") = snp_info,
		Rcpp::Named("targetID") = targetID,
		Rcpp::Named("genetype1") = genetype1,
		Rcpp::Named("genetype2") = genetype2,
		Rcpp::Named("expr_info") = expr_info,
		Rcpp::Named("out_param") = out_param);
		//Rcpp::Named("varcomp") = varcomp,
		//Rcpp::Named("h22") = h22);
		//Rcpp::Named("y_pred") = y_pred,
		//Rcpp::Named("y_se") = y_se);
		/*Rcpp::Named("tstat") = tstat,
		Rcpp::Named("alpha") = alpha,
		Rcpp::Named("sigma2beta") = sigma2betaout,
		Rcpp::Named("sigma2y") = sigma2yout,
		Rcpp::Named("sigma2z") = sigma2zout,
		Rcpp::Named("nsnpwgene") = nsnpwgene);*/
		//Rcpp::Named("idx") = idx);
		//Rcpp::Named("idxt3") = idxt3);

	return out;
}

