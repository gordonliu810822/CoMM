% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{CoMM_S4_testing}
\alias{CoMM_S4_testing}
\title{CoMM}
\usage{
CoMM_S4_testing(stringname1, stringname2, stringname3, stringname4,
  stringname5, px, lam)
}
\arguments{
\item{stringname1}{file for eQTL summary data file with the format described in the vignette.}

\item{stringname2}{file for GWAS summary data file with the format described in the vignette.}

\item{stringname3}{prefix for eQTL reference panel file with plink format (bim,bed).}

\item{stringname4}{prefix for the file complementing the information of SNPs in eQTL summary statistics with the format described in the vignette}

\item{stringname5}{prefix for GWAS reference panel file with plink format (bim,bed).}

\item{px}{logical value(1 or 0) for indicating whether to use paramter expansion to accerate algorithm}

\item{lam}{shirnkage intensity for reference panel}

}
\value{
List of model parameters
}
\description{
CoMM_S4_testing to dissecting genetic contributions to complex traits by leveraging regulatory information using both eQTL and GWAS summary statistics across whole genome 
}
\details{
\code{CoMM_S4_testing} fits the CoMM_S4 model across whole genome. It requires to provide both the eQTL and GWAS summary data, two reference panel data for eQTL and GWAS data respectively.
}
\examples{

file1 = "cis-eQTLs_full_chr1_1-100000.txt";
file2 = "NFBC_ph1_beta_se.txt";
file3 = "all_chr_1000G"
file4 = "cis-eQTLs_full_chr1";
file5 = "all_chr_1000G";
px = 1;
lam = 0.95;

fm = CoMM_S4_testing(file1,file2,file3, file4,file5, px, lam);

}
\author{
Yi Yang, \email{gmsyany@nus.edu.sg}
}
