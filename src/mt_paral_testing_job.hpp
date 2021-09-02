//  Created by Jin Liu, 27/04/2018.
//  Copyright © 2018年 Jin Liu. All rights reserved.
//

#ifndef mt_paral_testing_job_hpp
#define mt_paral_testing_job_hpp

#include <stdio.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <thread>
#include <mutex>
//#include "lmm_covar_pxem.hpp"
#include "lmm_covar_pxem_ptr.hpp"
//#include "AUDI_covar_pxem.hpp"
#include "CoMM_covar_pxem_ptr.hpp"
#include "plinkfun.hpp"

using namespace std;
using namespace arma;

class parGene_CoMM{
    public:
    int current_idx=0;

    // std::mutex mtx;

    //uvec bp_4use_r, chr_4use_r;
    //vec upper, lower, chr_expr, z;
    int N1, N2;
    long int P1, P2;
    uword Ngene_active;
    int bw;
    uword q;
	vec z, ind;
    
    char* X1;
    char* X2;
    mat expr_active, w1, w2, out_param, gene_info,snp_info;
    uvec idxinFile1, idxinFile2, idxin1;

    parGene_CoMM(char* X1, char* X2, const mat& expr_active, const vec& z, const mat& w1, const mat& w2,
		mat& out_param, const mat gene_info, const mat snp_info, const uword Ngene_active, const int bw, const int N1, const int N2, const long int P1, const long int P2, const uword q, const uvec& idxinFile1, const uvec& idxinFile2, const uvec& idxin1, const vec& ind){
        this -> X1 = X1;
        this -> X2 = X2;
        this -> expr_active = expr_active;
		this -> z = z;
        this -> w1 = w1;
		this -> w2 = w2;
        this -> out_param = out_param;
		this -> gene_info =gene_info;
		this -> snp_info = snp_info;
		this -> Ngene_active = Ngene_active;
        this -> bw = bw;
        this -> N1 = N1;
        this -> P1 = P1;
        this -> N2 = N2;
        this -> P2 = P2;
        this -> idxinFile1 = idxinFile1;
        this -> idxinFile2 = idxinFile2;
        this -> ind = ind;
        this -> idxin1 = idxin1;
		this -> q = q;
        //this -> idx_all = idx_all;
    }

    void loop_by_gene_CoMM(int i);
    void update_by_thread_CoMM(int thread_id);
    int  next_CoMM();

    /*mat X1;
    mat X2;
    mat X3;

    uword n;
    uword p;


    paraEx(const mat& X1, const mat& X2, mat& X3, uword n, uword p){

        this -> X1 = X1;
        this -> X2 = X2;
        this -> X3 = X3;
        this -> n = n;
        this -> p = p;

    }

    void loop_by_col(int i);
    void update_by_thread(int thread_id);
    int  next();*/

};



#endif /* mt_paral_job_hpp */

