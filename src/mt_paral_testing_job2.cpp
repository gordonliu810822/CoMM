#include "mt_paral_testing_job2.hpp"
#include "lmm_covar_pxem_ptr.hpp"
#include "SCLMM.hpp"



using namespace std;
using namespace arma;
using namespace Rcpp;



void parGene_SCLMM_IS::loop_by_gene_SCLMM(int g){

	uvec idx = find(snp_info.col(0) < gene_info(g, 2) + bw && snp_info.col(0) > gene_info(g, 1) - bw
		&& snp_info.col(1) == gene_info(g, 0));
	if (idx.n_elem == 0){
		cout << "Error: no matching SNPs for " << g << "-th gene ... " << endl;
	}

	// extract the SNPs within the region as X1 and X3
	mat X1tmp = X1.cols(idxinFile1(idx));
	X1tmp = X1tmp.rows(idxin1);
	mat X3tmp = X3.cols(idxinFile3(idx));

	uword idx1 = g;
	vec y = trans(expr_active.row(idx1));

	// remove the missing individuals from eQTL data
	uvec idx2 = find_finite(y);
	X1tmp = X1tmp.rows(idx2);
	y = y(idx2);

	// calculate cellular heritability
	rowvec sdX1tmp = stddev(X1tmp, 0, 0); // see manual
	mat X1tmp2 = X1tmp / repmat(sdX1tmp, X1tmp.n_rows, 1);

	// attention 
	uvec idx_finite = find_finite(X1tmp2.row(0));
	X1tmp2 = X1tmp2.cols(idx_finite);
	X1tmp = X1tmp.cols(idx_finite);
	X3tmp = X3tmp.cols(idx_finite);

	int Mg = X3tmp.n_cols;

	uvec idx_col1 = zeros<uvec>(1);
	uvec idx_col2 = ones<uvec>(1);
	vec hatmu2 = GWAS_sum(idxinFile2(idx), idx_col1);
	vec hats2 = GWAS_sum(idxinFile2(idx), idx_col2);

	// attention 
	hatmu2 = hatmu2.rows(idx_finite);
	hats2 = hats2.rows(idx_finite);

	mat R1 = eye(Mg, Mg);

	vec diagxtx3 = sum(X3tmp%X3tmp).t();
	for (int i = 0; i < Mg; i++){
		for (int j = 0; j < Mg; j++){
			if (i != j){
				R1(i, j) = sum(X3tmp.col(i) % X3tmp.col(j)) / sqrt(diagxtx3(i)*diagxtx3(j));
			}
		}
	}
	mat R;
	// positive definite sparse correlation matrix

	R = lam*R1 + (1 - lam)*eye(Mg, Mg);


	// set the parameter
	Options_SCLMM* lp_opt = new Options_SCLMM(100000, 0, 10, 1e-5, 0);
	Options_SCLMM* lp_opt1 = new Options_SCLMM(100000, 0, 10, 1e-5, 1);

	mat w1tmp = ones(X1tmp.n_rows, 1);

	double sigma2y, sigma2beta, loglik;
	vec beta0 = zeros<vec>(w1tmp.n_cols);
	int iterLmm, maxIter = 100000;
	mat Sigb = zeros<mat>(X1tmp.n_cols, X1tmp.n_cols);
	vec mub = zeros<vec>(X1tmp.n_cols);

	lmm_pxem_ptr(y, w1tmp, X1tmp2, maxIter, sigma2y, sigma2beta, beta0, loglik, iterLmm, Sigb, mub);

	out_param(g, 0) = sigma2beta*X1tmp.n_cols;
	out_param(g, 1) = sigma2y;
	out_param(g, 6) = 1 / (1 + sigma2y / (sigma2beta*X1tmp.n_cols));

	ObjSCLMM objHa = rcpparma_SCLMM_IS(X1tmp, y, w1tmp, hatmu2, hats2, R, lp_opt, 1);
	ObjSCLMM objH0 = rcpparma_SCLMM_IS(X1tmp, y, w1tmp, hatmu2, hats2, R, lp_opt1, 1);

	alphag(g) = objHa.alphag;
	stat(g) = 2 * (objHa.LRLB - objH0.LRLB);
	NumSNP(g) = Mg;

	out_param(g, 3) = 2 * (objHa.LRLB - objH0.LRLB);
	out_param(g, 4) = objHa.alphag;
	out_param(g, 5) = Mg;

	// reset
	X1tmp.reset();
	X1tmp2.reset();
	X3tmp.reset();
	w1tmp.reset();
	Sigb.reset();
	mub.reset();
	R1.reset();
	hatmu2.reset();
	hats2.reset();


	if ((g + 1) % 100 == 0 && (g + 1) != 0){
		cout << g + 1 << "-th Gene starts working ..." << endl;
	}

}

std::mutex _mtx22;
int parGene_SCLMM_IS::next_SCLMM(){
	std::lock_guard<std::mutex> lockGuard(_mtx22);
	if (current_idx >= Ngene){
		return -1;
	}
	current_idx++;
	return current_idx - 1;
}

void parGene_SCLMM_IS::update_by_thread_SCLMM(int thread_id){
	while (true){
		int idx = next_SCLMM();
		if (idx == -1){
			break;
		}
		loop_by_gene_SCLMM(idx);
	}
}


