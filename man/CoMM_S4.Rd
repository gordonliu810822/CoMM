% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{CoMM_S4}
\alias{CoMM_S4}
\title{CoMM}
\usage{
CoMM_S4(zscore, zscore2, R, R2, opts = NULL, px = 1)
}
\arguments{
\item{zscore}{zscore for eQTL summary statistics}

\item{zscore2}{zscore for GWAS summary statistics}

\item{R}{correlation matrix for LD information in eQTL data}

\item{R2}{correlation matrix for LD information in GWAS data}

\item{opts}{is a list consists of "max_iter", "dispF", "display_gap", "epsStopLogLik", "fix_alphag". 
        max_iter is an 
        positive integer which denotes the maximal iterations, default value 1e5. dispF is a logic value either 
        1 or 0, controlling whether or not 
        showing interations information, default value 1. display_gap is a positive integer which denotes the gap 
        showing interations information, default value 10. epsStopLogLik is a positive value which control 
        interation digits, default value 1e-5. constraintalpha is a logic value either 1 or 0, controlling whether
        alpha is constrained, default value 0. Noted that for separate analysis, fix_alphag does't exist in opts. 
        }

\item{px}{indicator for using PX-VBEM (default is 1). px = 1 for using PX-VBEM and VBEM otherwise.}
}
\value{
List of model parameters
}
\description{
CoMM_S4 to dissecting genetic contributions to complex traits by leveraging regulatory information using using both eQTL and GWAS summary data
}
\details{
\code{CoMM_S4} fits the TransCoMM_S2 model using both eQTL and GWAS summary data.
}
\examples{
set.seed(100)
L = 1;
M = 100;
rho =0.5;
n1 = 400;
n2 = 5000;
lam = 0.8;

X <- matrix(rnorm((n1+n2)*M),nrow=n1+n2,ncol=M);
beta_prop = 0.2;
b = matrix(0, M, 1); m = M * beta_prop;
b[sample(M,m)] = rnorm(m); h2y = 0.05;
y0 <- X\%*\%b + 6;
y  <- y0 + (as.vector(var(y0)*(1-h2y)/h2y))^0.5*rnorm(n1+n2);

h2 = 0.001;
y1 <- y[1:n1]
X1 <- X[1:n1,]
y2 <- y0[(n1+1):(n1+n2)]
X2 <- X[(n1+1):(n1+n2),]
alpha0 <- 3
alpha <- 0.3
sz2 <- var(y2*alpha) * ((1-h2)/h2)
z <- alpha0 + y2*alpha + rnorm(n2,0,sqrt(sz2))
y = y1;

mean.x1 = apply(X1,2,mean);
x1m = sweep(X1,2,mean.x1);
x1p = x1m;

mean.x2 = apply(X2,2,mean);
x2m = sweep(X2,2,mean.x2);
x2p = x2m;

w1 = matrix(rep(1,n1),ncol=1);

hatmu = matrix(0, M, 1)
hats = matrix(0, M, 1)

for (m in 1:M){
  fm = lm(y~1+x1p[,m]);
  hatmu[m] = summary(fm)$coefficients[2,1]
  hats[m] = summary(fm)$coefficients[2,2];
}
zscore = hatmu/hats;

hatmu2 = matrix(0, M, 1)
hats2 = matrix(0, M, 1)

for (m in 1:M){
  fm = lm(z~1+x2p[,m]);
  hatmu2[m] = summary(fm)$coefficients[2,1]
  hats2[m] = summary(fm)$coefficients[2,2];
}
zscore2 = hatmu2/hats2

x3p = x1p;              
sumx3p = apply(x3p*x3p, 2, sum)
R = matrix(0, M, M);
for (i1 in 1:M){
  for (j1 in 1:M){
    R[i1,j1] = t(x3p[,i1])\%*\%x3p[,j1]/sqrt(sumx3p[i1]*sumx3p[j1])                        
  }
}
R = R*lam + (1 - lam)*diag(M)  

x4p = x1p;              
sumx4p = apply(x4p*x4p, 2, sum)
R2 = matrix(0, M, M);
for (i1 in 1:M){
  for (j1 in 1:M){
    R2[i1,j1] = t(x4p[,i1])\%*\%x4p[,j1]/sqrt(sumx4p[i1]*sumx4p[j1])                        
  }
}
R2 = R2*lam + (1 - lam)*diag(M)  


opts = list(max_iter = 10000, dispF = 1, display_gap = 10, epsStopLogLik = 1e-5, fix_alphag = 0);
px = 1

fm = CoMM_S4(zscore, zscore2, R, R2, opts, px);

}
\author{
Yi Yang, \email{gmsyany@nus.edu.sg}
}
