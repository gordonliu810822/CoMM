#include "mt_paral_testing_job3.hpp"
#include "CoMM_S4.hpp"



using namespace std;
using namespace arma;
using namespace Rcpp;



void parGene_CoMM_2S2::loop_by_gene_CoMM_2S2(int g){

	//cout << "Info: checking point -3" << endl;
	int g_active = idx_active_gene(g);

	uvec idx31 = idx31_all(g_active);
	uvec idx32 = idxinFile2(idx_all(g_active));
	uvec idx33 = idxinFile3(idx_all(g_active));
	uvec idx35 = idxinFile5(idx_all(g_active));

	int Nsnp_gene = idx31.n_elem;

	//cout << "Info: checking point -2" << endl;

	uvec idx_col1 = zeros<uvec>(1);
	uvec idx_col2 = ones<uvec>(1);

	int Nsnp_gene_eQTL = gene_index(g_active, 1) - gene_index(g_active, 0) + 1;
	uvec idx41 = linspace<uvec>(gene_index(g_active, 0), gene_index(g_active, 1), Nsnp_gene_eQTL) - 1;
	uvec idx51 = idx41(idx31);

	//cout << "Info: checking point -1" << endl;

	vec hatmu1 = eqtls_sum(idx51, idx_col1);
	vec hats1 = eqtls_sum(idx51, idx_col2);
	vec hatmu2 = GWAS_sum(idx32, idx_col1);
	vec hats2 = GWAS_sum(idx32, idx_col2);
	vec zscore1 = hatmu1 / hats1;
	vec zscore2 = hatmu2 / hats2;

	//cout << "Info: checking point 0" << endl;

	cout << "Info: the number of snps for " << g << "-th gene is " << Nsnp_gene << endl;

	//First reference panel
	mat X_sub = X.cols(idx33);
	mat R_emp = eye(Nsnp_gene, Nsnp_gene);
	for (int i = 0; i < Nsnp_gene; i++){
		for (int j = 0; j < Nsnp_gene; j++){
			if (i != j){
				R_emp(i, j) = sum(X_sub.col(i) % X_sub.col(j)) / (norm(X_sub.col(i))*norm(X_sub.col(j)));
			}
		}
	}
	mat R = lam*R_emp + (1 - lam)*eye(R_emp.n_cols, R_emp.n_cols);

	//cout << "Info: checking point 1" << endl;

	//Second reference panel
	mat X2_sub = X2.cols(idx35);
	mat R2_emp = eye(Nsnp_gene, Nsnp_gene);
	for (int i = 0; i < Nsnp_gene; i++){
		for (int j = 0; j < Nsnp_gene; j++){
			if (i != j){
				R2_emp(i, j) = sum(X2_sub.col(i) % X2_sub.col(j)) / (norm(X2_sub.col(i))*norm(X2_sub.col(j)));
			}
		}
	}
	mat R2 = lam*R2_emp + (1 - lam)*eye(R2_emp.n_cols, R2_emp.n_cols);

	//cout << "Info: checking point 2" << endl;

	// set the parameter
	Options_CoMM_S4* lp_opt = new Options_CoMM_S4(10000, 0, 10, 1e-5, 0);
	Options_CoMM_S4* lp_opt1 = new Options_CoMM_S4(10000, 0, 10, 1e-5, 1);

	ObjCoMM_S4 objHa = rcpparma_CoMM_S4_VB(zscore1, zscore2, R, R2, lp_opt, px);
	ObjCoMM_S4 objH0 = rcpparma_CoMM_S4_VB(zscore1, zscore2, R, R2, lp_opt1, px);

	//cout << "Info: checking point 3" << endl;

	out_param(g, 0) = 2 * (objHa.LRLB - objH0.LRLB);
	out_param(g, 1) = objHa.alphag;
	out_param(g, 2) = Nsnp_gene;

	//cout << "Info: checking point 4" << endl;

	// reset
	X_sub.reset();
	X2_sub.reset();
	R_emp.reset();
	R2_emp.reset();
	hatmu1.reset();
	hats1.reset();
	hatmu2.reset();
	hats2.reset();



	if ((g + 1) % 100 == 0 && (g + 1) != 0){
		cout << g + 1 << "-th Gene starts working ..." << endl;
	}

}

std::mutex _mtx23;
int parGene_CoMM_2S2::next_CoMM_2S2(){
	std::lock_guard<std::mutex> lockGuard(_mtx23);
	if (current_idx >= Ngene_active){
		return -1;
	}
	current_idx++;
	return current_idx - 1;
}

void parGene_CoMM_2S2::update_by_thread_CoMM_2S2(int thread_id){
	while (true){
		int idx = next_CoMM_2S2();
		if (idx == -1){
			break;
		}
		loop_by_gene_CoMM_2S2(idx);
	}
}


