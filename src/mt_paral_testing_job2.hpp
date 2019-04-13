#ifndef mt_paral_testing_job_hpp
#define mt_paral_testing_job_hpp

#include <stdio.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <thread>
#include <mutex>
#include "SCLMM.hpp"
#include "lmm_covar_pxem_ptr.hpp"

using namespace std;
using namespace arma;
using namespace Rcpp;


class parGene_SCLMM_IS{
public:
	int current_idx = 0;
	int Ngene;

	mat gene_info, snp_info;
	int bw, N1, N3, P1, P3;
	uvec idxinFile1, idxinFile2, idxinFile3, idxin1;
	mat expr_active, GWAS_sum, X1, X3;
	vec alphag;
	vec stat;
	uvec NumSNP;
	double lam;
	string method;
	mat out_param;

	parGene_SCLMM_IS(mat& gene_info, mat& snp_info, int bw, int N1, int P1, int N3, int P3,
		uvec& idxinFile1, uvec& idxinFile2, uvec& idxinFile3, uvec& idxin1,
		mat& expr_active, mat& X1, mat& X3, mat& GWAS_sum,
		vec& alphag, vec& stat, uvec& NumSNP, int Ngene, double lam, string method, mat out_param){

		this->gene_info = gene_info;
		this->snp_info = snp_info;
		this->bw = bw;
		this->N1 = N1;
		this->P1 = P1;
		this->N3 = N3;
		this->P3 = P3;
		this->idxinFile1 = idxinFile1;
		this->idxinFile2 = idxinFile2;
		this->idxinFile3 = idxinFile3;
		this->idxin1 = idxin1;
		this->expr_active = expr_active;
		this->GWAS_sum = GWAS_sum;
		this->X1 = X1;
		this->X3 = X3;
		this->alphag = alphag;
		this->stat = stat;
		this->NumSNP = NumSNP;
		this->Ngene = Ngene;
		this->lam = lam;
		this->method = method;
		this->out_param = out_param;
	}

	void loop_by_gene_SCLMM(int g);
	void update_by_thread_SCLMM(int thread_id);
	int  next_SCLMM();

};


#endif 

