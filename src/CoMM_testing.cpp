#include <RcppArmadillo.h>
//#include <armadillo>
//#include <Rcpp.h>
//#include <omp.h>
#include <stdio.h>
#include <bitset>
#include <math.h>


#include "data_loader.hpp"
#include "mt_paral_testing_job.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;

//' @author Jin Liu, \email{jin.liu@duke-nus.edu.sg}
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
//' @param coreNum the number of core for parallelization. 
//'
//' @return List of model parameters
//'
//' @examples
//' file1 = "1000G.EUR.QC.1";
//' file2 = "NFBC_filter_mph10";
//' file3 = "Geuvadis_gene_expression_qn.txt";
//' file4 = "";
//' file5 = "pc5_NFBC_filter_mph10.txt";
//' whichPheno = 1;
//' bw = 500000;
//' coreNum = 24;
//'
//' fm = CoMM_testing_run_mt(file1,file2,file3, file4,file5, whichPheno, bw, coreNum);
//'
//' @details
//' \code{CoMM_testing_run_mt} fits the CoMM model. It requires to provide plink binary eQTL genotype file (bim, bed)
//' the GWAS plink binary file (bim, bed, fam), gene expression file for eQTL.
//' @export
// [[Rcpp::export]]
Rcpp::List CoMM_testing_run_mt(std::string stringname1, std::string stringname2, std::string stringname3, std::string stringname4,
                               std::string stringname5, int whCol, int bw,const int coreNum){ //int normalize_option = 1, int pred_option = 0){//, char* A21, char* A22){
    // normalize_option: 1. normalize each separately, 2. normalize both plink files together
    // match SNPs in file 1 and file 2 GWAS (common SNPs in x1 and x2 in columns)
    // plink file 1: stringname1; plink file 2: stringname2; expression file: stringname3
    // covariates file for file 1: stringname4; covariates file for file 2: stringname5
    // pred_option :0 (no calculation for prediction) 1 (calcuation for prediction)
    
    List tmp = dataLoader1(stringname1, stringname2, stringname3, whCol);
    //Mat<unsigned> X1 = tmp["X1"], X2 = tmp["X2"];
    vec z = tmp["y"];
    mat expr = tmp["expr_used"];
    CharacterVector rsname_4use_r = tmp["rsname_4use_r"];
    uvec chr_4use_r = tmp["chr_4use_r"];
    uvec bp_4use_r = tmp["bp_4use_r"];
    CharacterVector genetype1 = tmp["genetype1"], genetype2 = tmp["genetype2"], targetID = tmp["targetID"], indiv_4use = tmp["indiv_4use"];
    vec lower = tmp["lower"], upper = tmp["upper"], chr_expr = tmp["chr_expr"], ind = tmp["ind"];
    uvec idxin1 = tmp["idxin1"], idxin2 = tmp["idxin2"], idxinFile1 = tmp["idxinFile1"], idxinFile2 = tmp["idxinFile2"];
    
    // load the size of plink file
    string famfile1 = stringname1;
    famfile1 += ".fam";
    int N1 = getLineNum(famfile1);
    string bimfile1 = stringname1;
    bimfile1 += ".bim";
    long int P1 =  getLineNum(bimfile1);
    string famfile2 = stringname2;
    famfile2 += ".fam";
    int N2 = getLineNum(famfile2);
    string bimfile2 = stringname2;
    bimfile2 += ".bim";
    long int P2 =  getLineNum(bimfile2);
    
    long long size1 = (long long)N1 * (long long)P1;
    long long size2 = (long long)N2 * (long long)P2;
    
    char* X1 = new char[size1];
    char* X2 = new char[size2];
    char delimiter = '\t';
    
    clock_t t1 = clock();
    cout << "## Start loading genotype files 1, " ;
    readPlink2(stringname1,N1, P1, X1);
    cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
    
    t1 = clock();
    cout << "## Start loading genotype files 2, ";
    readPlink2(stringname2,N2, P2, X2);
    cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
    
    cout << "## Start loading covariates files ... " << endl;
    // load covariates file w2
    mat covar2;
    if (!stringname5.empty()){
        List tmp = getColNum_Header(stringname5, delimiter);
        int Ncovar = tmp["columns"];
        tmp = getCovarFile(stringname5, delimiter, Ncovar, N2);
        mat covartmp = tmp["covar"];
        mat w2one = ones<mat>(N2, 1);
        
        covar2 = join_rows(w2one, covartmp);
    }
    else {
        covar2 = ones<mat>(N2,1);
    }
    
	// load covariates file w1
	mat covar1;
	CharacterVector IID_w1;
	if (!stringname4.empty()){
		List tmp = getColNum_Header(stringname4, delimiter);

		int Ncovar = tmp["columns"];
		tmp = getCovarFile(stringname4, delimiter, Ncovar, idxin1.n_elem);
		mat covartmp = tmp["covar"];
		IID_w1 = tmp["IID"];
		IntegerVector idx_tmp = match(indiv_4use, IID_w1) -1;
		uvec idx_exprcov(as<uvec>(idx_tmp));

		covar1 = join_rows(ones<mat>(idxin1.n_elem, 1), covartmp.rows(idx_exprcov));
	}
	else {
		covar1 = ones<mat>(idxin1.n_elem,1);
	}
	cout << "## End loading files ... " << endl;
    
    mat w1 = covar1, w2 = covar2;
    
    uword q = w2.n_cols;
    
    //cout << "dim of w2 :" << w2.n_rows << "-by-" << q << endl;
    uword Ngene = lower.size();
    
    mat expr_info = zeros<mat>(lower.n_elem, 3);
    expr_info.col(0) = chr_expr;
    expr_info.col(1) = lower;
    expr_info.col(2) = upper;
    
    vec bp_4use_r_vec = conv_to<vec>::from(bp_4use_r);
    vec chr_4use_r_vec = conv_to<vec>::from(chr_4use_r);
    
    mat snp_info = zeros<mat>(chr_4use_r.n_elem, 2);
    snp_info.col(0) = bp_4use_r_vec;
    snp_info.col(1) = chr_4use_r_vec;
    
    uvec idx;
    
    //uvec idx_used ;
    t1 = clock();
    
    cout << "### Precompute the active genes ... ";
    //uvec idx_all = zeros<uvec>(0);
    uvec idx_active_gene = zeros<uvec>(0);
    uvec g_tmp(1);
    for (uword g = 0; g < Ngene; g++){
        
        idx = find(bp_4use_r_vec < upper(g) + bw && bp_4use_r_vec > lower(g) - bw
                   && chr_4use_r_vec == chr_expr(g));
        if (idx.is_empty() == false){
            if (idx.n_elem > 1){
                //idx_all = join_cols(idx_all, idx);
                g_tmp(0) = g;
                idx_active_gene = join_cols(idx_active_gene,g_tmp);
            }
        }
    }
    cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
    //cout << "dim of out_param0" << out_param0.n_rows <<"-by-" << out_param0.n_cols << endl;
    
    // combine gene info into mat
    uword Ngene_active = idx_active_gene.n_elem;
    mat gene_info = zeros<mat>(Ngene_active,3);
    gene_info = expr_info.rows(idx_active_gene);
    mat expr_active = expr.rows(idx_active_gene);
    
    CharacterVector gene_type1 = genetype1[as<Rcpp::IntegerVector>(wrap(idx_active_gene))];
    CharacterVector gene_type2 = genetype2[as<Rcpp::IntegerVector>(wrap(idx_active_gene))];
    
    mat out_param0 = -99*ones<mat>(7, Ngene_active);
    //out_param0 = -99*out_param0;
    
    cout << "### Start running CoMM on active genes ... " << endl;;
    t1 = clock();
    //set parallele structure object
    parGene_CoMM parObj(X1,X2,expr_active,z,w1,w2,out_param0,gene_info,snp_info,Ngene_active,bw,N1,N2,P1,P2,q,idxinFile1,idxinFile2,idxin1,ind);
    
    //set parallel computation
    const int n_thread = coreNum;
    std::vector<std::thread> threads(n_thread);
    for(int i_thread = 0; i_thread < n_thread; i_thread++){
        threads[i_thread] = std::thread(&parGene_CoMM::update_by_thread_CoMM, &parObj, i_thread);
    }
    
    for(int i = 0; i < n_thread; i++){
        threads[i].join();
    }
    
    delete[] X1;
    delete[] X2;
    
    cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
    
    mat out_param = -99*ones<mat>(Ngene,7);
    out_param.rows(idx_active_gene) = trans(parObj.out_param);
    
    List out = List::create(
                            Rcpp::Named("z") = z,
                            Rcpp::Named("covar1") = w1,
                            Rcpp::Named("covar2") = w2,
                            Rcpp::Named("expr") = expr,
                            Rcpp::Named("expr_active") = expr_active,
                            Rcpp::Named("rsname_4use_r") = rsname_4use_r,
                            //Rcpp::Named("chr_4use_r") = chr_4use_r,
                            //Rcpp::Named("bp_4use_r") = bp_4use_r,
                            Rcpp::Named("snp_info_4use_r") = snp_info,
                            Rcpp::Named("targetID") = targetID,
                            Rcpp::Named("genetype1") = genetype1,
                            Rcpp::Named("genetype2") = genetype2,
                            Rcpp::Named("gene_type1") = gene_type1,
                            Rcpp::Named("gene_type2") = gene_type2,
                            Rcpp::Named("gene_info") = gene_info,
                            Rcpp::Named("expr_info") = expr_info,
                            Rcpp::Named("out_param") = out_param,
                            Rcpp::Named("out_param0") = parObj.out_param);
    
    return out;
}

//' @author Jin Liu, \email{jin.liu@duke-nus.edu.sg}
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
//' file1 = "1000G.EUR.QC.1";
//' file2 = "NFBC_filter_mph10";
//' file3 = "Geuvadis_gene_expression_qn.txt";
//' file4 = "";
//' file5 = "pc5_NFBC_filter_mph10.txt";
//' whichPheno = 1;
//' bw = 500000;
//'
//' fm = CoMM_testing_run(file1,file2,file3, file4,file5, whichPheno, bw);
//'
//' @details
//' \code{CoMM_testing_run} fits the CoMM model. It requires to provide plink binary eQTL genotype file (bim, bed)
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
    
    List tmp = dataLoader1(stringname1, stringname2, stringname3, whCol);
    vec z = tmp["y"];
    mat expr = tmp["expr_used"];
    CharacterVector rsname_4use_r = tmp["rsname_4use_r"];
    uvec chr_4use_r = tmp["chr_4use_r"];
    uvec bp_4use_r = tmp["bp_4use_r"];
    CharacterVector genetype1 = tmp["genetype1"], genetype2 = tmp["genetype2"], targetID = tmp["targetID"], indiv_4use = tmp["indiv_4use"];
    vec lower = tmp["lower"], upper = tmp["upper"], chr_expr = tmp["chr_expr"], ind = tmp["ind"];
    uvec idxin1 = tmp["idxin1"], idxin2 = tmp["idxin2"], idxinFile1 = tmp["idxinFile1"], idxinFile2 = tmp["idxinFile2"];
    
    // load the size of plink file
    string famfile1 = stringname1;
    famfile1 += ".fam";
    int N1 = getLineNum(famfile1);
    string bimfile1 = stringname1;
    bimfile1 += ".bim";
    long int P1 =  getLineNum(bimfile1);
    string famfile2 = stringname2;
    famfile2 += ".fam";
    int N2 = getLineNum(famfile2);
    string bimfile2 = stringname2;
    bimfile2 += ".bim";
    long int P2 =  getLineNum(bimfile2);
    
    long long size1 = (long long)N1 * (long long)P1;
    long long size2 = (long long)N2 * (long long)P2;
    
    char* X1 = new char[size1];
    char* X2 = new char[size2];
    char delimiter = '\t';
    
    clock_t t1 = clock();
    cout << "## Start loading genotype files 1, " ;
    readPlink2(stringname1,N1, P1, X1);
    cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
    
    t1 = clock();
    cout << "## Start loading genotype files 2, ";
    readPlink2(stringname2,N2, P2, X2);
    cout << ", Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
    
    cout << "## Start loading covariates files ... " << endl;
    // load covariates file w2
    mat covar2;
    if (!stringname5.empty()){
        List tmp = getColNum_Header(stringname5, delimiter);
        int Ncovar = tmp["columns"];
        tmp = getCovarFile(stringname5, delimiter, Ncovar, N2);
        mat covartmp = tmp["covar"];
        mat w2one = ones<mat>(N2, 1);
        
        covar2 = join_rows(w2one, covartmp);
    }
    else {
        covar2 = ones<mat>(N2,1);
    }
    
	// load covariates file w1
	mat covar1;
	CharacterVector IID_w1;
	if (!stringname4.empty()){
		List tmp = getColNum_Header(stringname4, delimiter);

		int Ncovar = tmp["columns"];
		tmp = getCovarFile(stringname4, delimiter, Ncovar, idxin1.n_elem);
		mat covartmp = tmp["covar"];
		IID_w1 = tmp["IID"];
		IntegerVector idx_tmp = match(indiv_4use, IID_w1) -1;
		uvec idx_exprcov(as<uvec>(idx_tmp));

		covar1 = join_rows(ones<mat>(idxin1.n_elem, 1), covartmp.rows(idx_exprcov));
	}
	else {
		covar1 = ones<mat>(idxin1.n_elem,1);
	}
	cout << "## End loading files ... " << endl;
    
    mat w1 = covar1, w2 = covar2;
    
    uword Ngene = lower.size();
    umat snp_info = zeros<umat>(chr_4use_r.n_elem, 2);
    snp_info.col(0) = chr_4use_r;
    snp_info.col(1) = bp_4use_r;
    
    mat expr_info = zeros<mat>(lower.n_elem, 3);
    expr_info.col(0) = chr_expr;
    expr_info.col(1) = lower;
    expr_info.col(2) = upper;
    
    uvec idx;
    
    vec bp_4use_r_vec = conv_to<vec>::from(bp_4use_r);
    vec chr_4use_r_vec = conv_to<vec>::from(chr_4use_r);
    
    mat out_param = -99*ones<mat>(Ngene,7);
    
    uvec idx_all = zeros<uvec>(0);
    uvec idx_active_gene = zeros<uvec>(0);
    uvec g_tmp(1);
    int maxIter = 1000;
    double epsStopLogLik = 1e-5;
    int constr;
    
    // conduct testing
    
    t1 = clock();
    for (uword g = 0; g < Ngene; g++){
        
        idx = find(bp_4use_r_vec < upper(g) + bw && bp_4use_r_vec > lower(g) - bw
                   && chr_4use_r_vec == chr_expr(g));
        
        
        out_param(g, 5) = idx.n_elem;
        
        
        if (idx.is_empty() == false){
            if (idx.n_elem > 1){
                g_tmp(0) = g;
                idx_active_gene = join_cols(idx_active_gene,g_tmp);
                
                
                if ( idx_active_gene.n_elem % 100 == 0 && idx_active_gene.n_elem != 0){
                    cout << idx_active_gene.n_elem << "-th Gene starts working ..." ;
                    cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
                }
                
                idx_all = join_cols(idx_all, idx);
                
                double* sub_matrix_double1 = new double[idx.n_elem * N1];
                double* sub_matrix_double2 = new double[idx.n_elem * N2];
                
                mat X1tmp = getSubMat(X1, N1 , P1, idxinFile1(idx)-1, sub_matrix_double1);
                X1tmp = X1tmp.rows(idxin1);
                X1tmp.replace(3, 0);
                mat X2tmp = getSubMat(X2, N2 , P2, idxinFile2(idx)-1, sub_matrix_double2);
                X2tmp.replace(3, 0);
                
                
                vec ind_idx = ind(idx);
                uvec ind_idx_1 = find(ind_idx == -1);
                
                X2tmp.cols(ind_idx_1) = 2 - X2tmp.cols(ind_idx_1);
                
                uvec idx1(1); idx1(0) = g;
                vec y = trans(expr.rows(idx1));
                y.replace(datum::nan, 9999);
                
                uvec idx2 = find(y != 9999);
                X1tmp = ((&X1tmp) -> rows(idx2));
                
                
                mat w1tmp = ((&w1) -> rows(idx2));
                
                y = y(idx2);
            
                rowvec meanX1tmp = mean(X1tmp, 0);
                rowvec sdX1tmp = stddev(X1tmp, 0, 0); // see manual
                X1tmp = (X1tmp - repmat(meanX1tmp, X1tmp.n_rows, 1))/ repmat(sdX1tmp, X1tmp.n_rows, 1) / sqrt(X1tmp.n_cols);
                
                rowvec meanX2tmp = mean(X2tmp, 0);
                rowvec sdX2tmp = stddev(X2tmp, 0, 0); // see manual
                X2tmp = (X2tmp - repmat(meanX2tmp, X2tmp.n_rows, 1))/ repmat(sdX2tmp, X2tmp.n_rows, 1) / sqrt(X2tmp.n_cols);
                
                rowvec X1row = X1tmp.row(0);
                X1tmp = X1tmp.cols(find_finite(X1row));
                
                X2tmp = X2tmp.cols(find_finite(X1row));
                
                // initialize by linear mixed model for sigma2y and sigma2beta
                // use this initialization to calculate MoM estimates of h2
                double sigma2y, sigma2beta, loglik;
                vec beta0 =zeros<vec>(w1tmp.n_cols);
                int iter;
                mat Sigb = zeros<mat>(X1tmp.n_cols,X1tmp.n_cols);
                vec mub  = zeros<vec>(X1tmp.n_cols);
                
                lmm_pxem_ptr(y, w1tmp, X1tmp, maxIter,sigma2y,sigma2beta,beta0,loglik,iter,Sigb,mub);
                
                out_param(g, 6) = 1/(1+sigma2y/sigma2beta); // gene-wise h2
                
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
                
                //remove local variable by reset
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
                
                out_param(g, 0) = sigma2betaHa;
                out_param(g, 1) = sigma2yHa;
                out_param(g, 2) = sigma2zHa;
                out_param(g, 3) = 2 * (loglikHa - loglikHo);
                out_param(g, 4) = alphaHa;
            }
            else{
            }
        }
        
        else{
        }
        
    }
    delete[] X1;
    delete[] X2;
    
    cout << "Model fitting is done in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
    
    // combine gene info into mat
    uword Ngene_active = idx_active_gene.n_elem;
    mat gene_info = zeros<mat>(Ngene_active,3);
    gene_info = expr_info.rows(idx_active_gene);
    mat out_param0 = out_param.rows(idx_active_gene);
    
    CharacterVector gene_type1 = genetype1[as<Rcpp::IntegerVector>(wrap(idx_active_gene))];
    CharacterVector gene_type2 = genetype2[as<Rcpp::IntegerVector>(wrap(idx_active_gene))];
    
    
    
    
    List out = List::create(Rcpp::Named("idx_all") = idx_all,
                            Rcpp::Named("z") = z,
                            Rcpp::Named("covar1") = w1,
                            Rcpp::Named("covar2") = w2,
                            Rcpp::Named("expr") = expr,
                            Rcpp::Named("rsname_4use_r") = rsname_4use_r,
                            Rcpp::Named("snp_info_4use_r") = snp_info,
                            Rcpp::Named("targetID") = targetID,
                            Rcpp::Named("genetype1") = genetype1,
                            Rcpp::Named("genetype2") = genetype2,
                            Rcpp::Named("gene_type1") = gene_type1,
                            Rcpp::Named("gene_type2") = gene_type2,
                            Rcpp::Named("expr_info") = expr_info,
                            Rcpp::Named("gene_info") = gene_info,
                            Rcpp::Named("out_param") = out_param,
                            Rcpp::Named("out_param0") = out_param0);
    
    return out;
}


