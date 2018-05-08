//  Created by Jin Liu, 27/04/2018.
//  Copyright © 2018年 Jin Liu. All rights reserved.
//


#include "mt_paral_testing_job.hpp"

using namespace std;
using namespace arma;

void parGene_CoMM::loop_by_gene_CoMM(int i){

    vec out_param_par = vec(out_param.colptr(i),7, false);

    uvec idx = find(snp_info.col(0) < gene_info(i,2) + bw && snp_info.col(0) > gene_info(i,1) - bw
        && snp_info.col(1) == gene_info(i,0));
    if (idx.n_elem == 0 ){
        cout << "Error: no matching SNPs for " << i << "-th gene ... " << endl;
    }
    out_param_par(5) = idx.n_elem;

    // extract the SNPs within the region as X1 and X2
	//mat X1tmp = conv_to<mat>::from(X1.cols(idx));//conv_to<mat>::from((&X1)->cols(idx));
	//mat X2tmp = conv_to<mat>::from(X2.cols(idx));//conv_to<mat>::from((&X2)->cols(idx));
    double* sub_matrix_double1 = new double[idx.n_elem * N1];
    double* sub_matrix_double2 = new double[idx.n_elem * N2];

    mat X1tmp = getSubMat(X1, N1 , P1, idxinFile1(idx)-1, sub_matrix_double1);
    X1tmp = X1tmp.rows(idxin1);
    X1tmp.replace(3, 0);
    mat X2tmp = getSubMat(X2, N2 , P2, idxinFile2(idx)-1, sub_matrix_double2);
    // X2tmp = X2tmp.rows(idxin2);
    X2tmp.replace(3, 0);

    //cout << "Dim of X2tmp: " << X2tmp.n_rows <<"*" <<X2tmp.n_cols << endl;

    vec ind_idx = ind(idx);
    uvec ind_idx_1 = find(ind_idx == -1);

    X2tmp.cols(ind_idx_1) = 2 - X2tmp.cols(ind_idx_1);
	//cout << i << "-th gene: sum(X1tmp) " << sum(sum(X1tmp)) << "; sum(X2tmp) " << sum(sum(X2tmp)) << endl;
	//extract i-th row from expr
    uword idx1 = i;
	vec y = trans(expr_active.row(idx1));

    // remove the missing individuals from eQTL data
    uvec idx2 = find_finite(y);
	X1tmp = X1tmp.rows(idx2);//((&X1tmp)->rows(idx2));
	y = y(idx2);

	mat w1tmp = w1.rows(idx2);//((&w1)->rows(idx2));

	rowvec meanX1tmp = mean(X1tmp, 0);
	rowvec sdX1tmp = stddev(X1tmp, 0, 0); // see manual
	rowvec meanX2tmp = mean(X2tmp, 0);
	rowvec sdX2tmp = stddev(X2tmp, 0, 0); // see manual

    // normalize X1 and X2
	X1tmp = (X1tmp - repmat(meanX1tmp, X1tmp.n_rows, 1)) / repmat(sdX1tmp, X1tmp.n_rows, 1)/ sqrt(X1tmp.n_cols);;

	X2tmp = (X2tmp - repmat(meanX2tmp, X2tmp.n_rows, 1))/ repmat(sdX2tmp, X2tmp.n_rows, 1) / sqrt(X2tmp.n_cols);

	rowvec X1row = X1tmp.row(0);

    uvec idx3 = find_finite(X1row);
	X1tmp = X1tmp.cols(idx3);

	X2tmp = X2tmp.cols(idx3);

    if (idx3.n_elem == 0){
        cout << "Error: the number of SNPs are 0 for " << i << "-th gene ... " << endl;
    }

	double sigma2y, sigma2beta, loglik;
	vec beta0 =zeros<vec>(w1tmp.n_cols);
	int iterLmm, maxIter=1000;
	mat Sigb = zeros<mat>(X1tmp.n_cols,X1tmp.n_cols);
	vec mub  = zeros<vec>(X1tmp.n_cols);
	lmm_pxem_ptr(y, w1tmp, X1tmp, maxIter,sigma2y,sigma2beta,beta0,loglik,iterLmm,Sigb,mub);

	out_param_par(6) = 1/(1+sigma2y/sigma2beta);

	double epsStopLogLik = 1e-5;
	int constr;

	double sigma2betaHa = sigma2beta, sigma2betaHo = sigma2beta, sigma2yHa = sigma2y, sigma2yHo = sigma2y;
	double sigma2zHa, sigma2zHo, gamHa, gamHo, loglikHa, loglikHo, alphaHa, alphaHo;
	int iterHa, iterHo;
	vec alpha0Ha = zeros<vec>(w2.n_cols), alpha0Ho = zeros<vec>(w2.n_cols);
	vec beta0Ha = beta0, beta0Ho = beta0;

	//Ha: run
	constr = 0;
	CoMM_covar_pxem_ptr(y, z, X1tmp, X2tmp, w1tmp, w2, sigma2betaHa, sigma2yHa, sigma2zHa, alpha0Ha, alphaHa, gamHa, beta0Ha, loglikHa,
		iterHa, constr, epsStopLogLik, maxIter);
	//H0: run
	constr = 1;
	CoMM_covar_pxem_ptr(y, z, X1tmp, X2tmp, w1tmp, w2, sigma2betaHo, sigma2yHo, sigma2zHo, alpha0Ho, alphaHo, gamHo, beta0Ho, loglikHo,
		iterHo, constr, epsStopLogLik, maxIter);

    //delete[] X1tmp, X2tmp, w1tmp;
    X1tmp.reset();
    X2tmp.reset();
    w1tmp.reset();
    Sigb.reset();
    mub.reset();
    meanX1tmp.reset();
    meanX2tmp.reset();
    sdX1tmp.reset();
    sdX2tmp.reset();
    X1row.reset();
    ind_idx.reset();
    ind_idx_1.reset();
    w1tmp.reset();
    delete[] sub_matrix_double1;
    delete[] sub_matrix_double2;

	out_param_par(0) = sigma2betaHa;
	out_param_par(1) = sigma2yHa;
	out_param_par(2) = sigma2zHa;
	out_param_par(3) = 2 * (loglikHa - loglikHo); //2 * ( loglikHa(loglikHa.n_elem-1) - loglikHo(loglikHo.n_elem-1));//
	out_param_par(4) = alphaHa;

    if ( (i+1) % 100 == 0 && (i+1) != 0){
        cout << i+1 << "-th Gene starts working ..." << endl;
    }
	//}
	//cout << i << "-th gene: sigma2beta: "<< sigma2betaHa << "; sigma2yHa: " << sigma2yHa << "; sigma2z: " << sigma2zHa<< "; lrt: " << out_param_par(3) << endl;
	// }
    //}

}

std::mutex _mtx2;
int parGene_CoMM::next_CoMM(){
    std::lock_guard<std::mutex> lockGuard(_mtx2);
    if(current_idx >= (int)Ngene_active){
        return -1;
    }
    current_idx++;
    return current_idx-1;
}

void parGene_CoMM::update_by_thread_CoMM(int thread_id){
    while(true){
        int idx = next_CoMM();
        if(idx == -1){
            break;
        }
        // cout << idx << endl;
        loop_by_gene_CoMM(idx);
    }
}

