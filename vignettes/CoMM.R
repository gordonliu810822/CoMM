## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(devtools)
install_github("gordonliu810822/CoMM")

## ------------------------------------------------------------------------
library("CoMM")

## ------------------------------------------------------------------------
L = 1; M = 100; rho =0.5
n1 = 350; n2 = 5000;
maf = runif(M,0.05,0.5)
X = genRawGeno(maf, L, M, rho, n1 + n2);

## ------------------------------------------------------------------------
beta_prop = 0.2;
b = numeric(M);
m = M * beta_prop;
b[sample(M,m)] = rnorm(m);

## ------------------------------------------------------------------------
h2y = 0.05;
b0 = 6;
y0 <- X%*%b + b0;
y  <- y0 + (as.vector(var(y0)*(1-h2y)/h2y))^0.5*rnorm(n1+n2);

## ------------------------------------------------------------------------
h2 = 0.001;
y1 <- y[1:n1]
X1 <- X[1:n1,]
y2 <- y0[(n1+1):(n1+n2)]
X2 <- X[(n1+1):(n1+n2),]
alpha0 <- 3  
alpha <- 0.3
sz2 <- var(y2*alpha) * ((1-h2)/h2)
z <- alpha0 + y2*alpha + rnorm(n2,0,sqrt(sz2))

## ------------------------------------------------------------------------
y = y1;
mean.x1 = apply(X1,2,mean);
x1m = sweep(X1,2,mean.x1);
std.x1 = apply(x1m,2,sd)
x1p = sweep(x1m,2,std.x1,"/");
x1p = x1p/sqrt(dim(x1p)[2])

mean.x2 = apply(X2,2,mean);
x2m = sweep(X2,2,mean.x2);
std.x2 = apply(x2m,2,sd)
x2p = sweep(x2m,2,std.x2,"/");
x2p = x2p/sqrt(dim(x2p)[2])

w2 = matrix(rep(1,n2),ncol=1);
w1 = matrix(rep(1,n1),ncol=1);

## ------------------------------------------------------------------------
fm0 = lmm_pxem(y, w1,x1p, 100)
sigma2beta =fm0$sigma2beta;
sigma2y =fm0$sigma2y;
beta0 = fm0$beta0;

## ------------------------------------------------------------------------
fmHa = CoMM_covar_pxem(y, z, x1p, x2p, w1, w2, sigma2beta, sigma2y, beta0, 0, 1e-5, 1000);
fmH0 = CoMM_covar_pxem(y, z, x1p, x2p, w1, w2, sigma2beta, sigma2y, beta0, 1, 1e-5, 1000);
loglikHa = max(fmHa$loglik,na.rm=T)
loglikH0 = max(fmH0$loglik,na.rm=T)
tstat = 2 * (loglikHa - loglikH0);
pval = pchisq(tstat,1,lower.tail=F)
alpha_hat = fmHa$alpha 

## ------------------------------------------------------------------------
file1 = "1000G.EUR.QC.1";
file2 = "NFBC_filter_mph10";
file3 = "Geuvadis_gene_expression_qn.txt";
file4 = "";
file5 = "pc5_NFBC_filter_mph10.txt";
whichPheno = 1;
bw = 500000;

## ----table2, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "  # simple table creation here
| Tables        | Are           | Cool  |
|---------------|:-------------:|------:|
| col 3 is      | right-aligned | $1600 |
| col 2 is      | centered      |   $12 |
| zebra stripes | are neat      |    $1 |
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

