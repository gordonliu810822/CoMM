#ifndef mt_paral_testing_job3_hpp
#define mt_paral_testing_job3_hpp

#include <stdio.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <thread>
#include <mutex>
#include "CoMM_S4.hpp"

using namespace std;
using namespace arma;
using namespace Rcpp;


class parGene_CoMM_2S2{
public:
	int current_idx = 0;
	int Ngene_active;

	uvec idxinFile5, idxinFile2, idxinFile3, idxin1;
	mat GWAS_sum, eqtls_sum, X, X2;
	uvec NumSNP;
	mat out_param;
	umat gene_index;
	field<uvec> idx31_all, idx_all;
	bool px;
	uvec idx_active_gene;
	double lam;

	parGene_CoMM_2S2(umat& gene_index, field<uvec>& idx31_all, field<uvec>& idx_all, uvec& idxinFile2, uvec& idxinFile3, uvec& idxinFile5, mat& GWAS_sum, mat& eqtls_sum,
		mat& X, mat& X2, bool px, mat out_param, int Ngene_active, uvec idx_active_gene, double lam){

		this->gene_index = gene_index;
		this->idx31_all = idx31_all;
		this->idx_all = idx_all;
		this->idxinFile2 = idxinFile2;
		this->idxinFile3 = idxinFile3;
		this->idxinFile5 = idxinFile5;
		this->GWAS_sum = GWAS_sum;
		this->eqtls_sum = eqtls_sum;
		this->X = X;
		this->X2 = X2;
		this->px = px;
		this->out_param = out_param;
		this->Ngene_active = Ngene_active;
		this->idx_active_gene = idx_active_gene;
		this->lam = lam;
	}

	void loop_by_gene_CoMM_2S2(int g);
	void update_by_thread_CoMM_2S2(int thread_id);
	int  next_CoMM_2S2();

};


#endif 

