## file1 corresponds to the individual-level transcriptome data
## file2 corresponds to the gene expression file
## file3 corresponds to the individual-level GWAS data

## file4 corresponds to summary-level transcriptome data
## file5 corresponds to complementary summary-level transcriptome data
## file6 corresponds to summary-level GWAS data

## Cov1 corresponds to corvariates for transcriptome
## Cov2 corresponds to corvariates for GWAS

## R1 corresponds to reference panel for transcriptome
## R2 corresponds to reference panel for GWAS

#' CoMM.
#' 
#' @description
#' CoMM implement the model CoMM, CoMM-S2 and CoMM-S4. 
#'
#' @details CoMM implement the model CoMM, CoMM-S2 and CoMM-S4. 
#' @param Model an character specifying which model to run. "Individual" indicates the model using individual eQTL and GWAS data.
#' "Summary" indicates the model using individual eQTL and summary GWAS data. "Trans" indicates the model using both eQTL and GWAS data
#' with different reference panel for eQTL and GWAS data. "Shared" indicates the model using both eQTL and GWAS data with the shared reference panel.
#' @param file1 is the name of transcriptome genotype file.
#' @param file2 is the name of transcriptome gene expression file.
#' @param file3 is the name of GWAS plink file.
#' @param file4 is the name of the summary-level transcriptome data.
#' @param file5 is the name of complementary summary-level transcriptome data.
#' @param file6 is the name of summary-level GWAS data.
#' @param Cov1 is the name of corvariates for transcriptome.
#' @param Cov2 is the name of corvariates for GWAS.
#' @param R1 is the name of reference panel for transcriptome.
#' @param R2 is the name of reference panel for GWAS.
#' @param whCol is an integer specifying which column of phenotype. The default is 1.
#' @param bw is an integer specifying the number of downstream and upstream SNPs that are considered as cis-SNP within a gene. The default is 500000.
#' @param lam is a positive number specifying the shrinkage intensity for reference panel. The default is 0.95.
#' @param paral is a logical value specifying whether to run the model in parallel. The default is FALSE.
#' @param coreNum is an integer specifying the number of cores for parallel. The default is 1.
#' @return a list, We briefly explain the output of the CoMM. 
#'  
#' @export
CoMM <- function(Model = "Individual", file1 = NULL, file2 = NULL, file3 = NULL, file4 = NULL, file5 = NULL, file6 = NULL, Cov1 = NULL, Cov2 = NULL, R1 = NULL, R2 = NULL, whCol = 1, bw = 500000, lam = 0.95, paral = FALSE, coreNum = 1){
  
  ## Model Individual
  if (Model == "Individual"){
    if (is.null(file1)){
      stop("Require transcriptome genotype file: file1")
    }
    if (is.null(file2)){
      stop("Require transcriptome gene expression file: file2")
    }
    if (is.null(file3)){
      stop("Require GWAS plink file: file3")
    }
    if (paral == TRUE){
      if (is.null(Cov1)){
        Cov1 = ""
      }
      if (is.null(Cov2)){
        Cov2 = ""
      }
      .Call('_CoMM_CoMM_testing_run_mt', PACKAGE = 'CoMM', file1, file3, file2, Cov1, Cov2, whCol, bw, coreNum)
    }else{
      if (is.null(Cov1)){
        Cov1 = ""
      }
      if (is.null(Cov2)){
        Cov2 = ""
      }
      .Call('_CoMM_CoMM_testing_run', PACKAGE = 'CoMM', file1, file3, file2, Cov1, Cov2, whCol, bw)
    }
    
  }else if (Model == "Summary"){
    if (is.null(file1)){
      stop("Require transcriptome genotype file: file1")
    }
    if (is.null(file2)){
      stop("Require transcriptome gene expression file: file2")
    }
    if (is.null(file6)){
      stop("Require summary-level GWAS file: file6")
    }
    if (is.null(R2)){
      stop("Require reference panel for GWAS: R2")
    }
    
    if (paral == TRUE){
      if (is.null(Cov1)){
        Cov1 = ""
      }
     .Call('_CoMM_CoMM_S2_paral_testing', PACKAGE = 'CoMM', file1, file6, R2, file2, Cov1, bw, lam, coreNum)
    }else{
      if (is.null(Cov1)){
        Cov1 = ""
      }
      .Call('_CoMM_CoMM_S2_testing', PACKAGE = 'CoMM', file1, file6, R2, file2, Cov1, bw, lam, 0)
    }
    
  }else if (Model == "Trans"){
    if (is.null(file4)){
      stop("Require summmary-level eQTL file: file4")
    }
    if (is.null(file5)){
      stop("Require complementary summary-level eQTL file: file5")
    }
    if (is.null(file6)){
      stop("Require summary-level GWAS file: file6")
    }
    if (is.null(R1)){
      stop("Require reference panel for eQTL: R1")
    }
    if (is.null(R2)){
      stop("Require reference panel for GWAS: R2")
    }
    if (paral == TRUE){
      .Call('_CoMM_CoMM_S4_testing_mt', PACKAGE = 'CoMM', file4, file6, R1, file5, R2, 1, lam, coreNum)
    }else{
      .Call('_CoMM_CoMM_S4_testing', PACKAGE = 'CoMM', file4, file6, R1, file5, R2, 1, lam)
    }
    
  }else if (Model == "Shared"){
    if (is.null(file4)){
      stop("Require summmary-level eQTL file: file4")
    }
    if (is.null(file5)){
      stop("Require complementary summary-level eQTL file: file5")
    }
    if (is.null(file6)){
      stop("Require summary-level GWAS file: file6")
    }
    if (is.null(R1)){
      stop("Require reference panel both for transcriptome and GWAS: R1")
    }
    if (paral == TRUE){
      .Call('_CoMM_CoMM_S4_testing_mt', PACKAGE = 'CoMM', file4, file6, R1, file5, R1, 1, lam, coreNum)
    }else{
      .Call('_CoMM_CoMM_S4_testing', PACKAGE = 'CoMM', file4, file6, R1, file5, R1, 1, lam)
    }
  }else{
    print("Wrong Model name, please specify the right Model name: Individual, Summary, Trans or Shared")
  }

}
  


