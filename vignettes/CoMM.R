## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library("CoMM")

## ------------------------------------------------------------------------
library(mvtnorm)
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

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(ge[,1:8])

## ---- fig.width=7, fig.height=4, echo=FALSE------------------------------
dat_rej = dat[[3]];
dat_rej$h2z=paste("",dat_rej$h2,sep="")
dat_rej$Power = dat_rej$rej_prop
dat_rej$Sparsity = dat_rej$beta_prop
dat_rej$sd_rej = as.numeric(as.character(dat_rej$sd_rej))
dat_rej = dat_rej[dat_rej$Method!="2-stage:AUDI",]
library(plyr)
dat_rej$Method=revalue(dat_rej$Method, c("AUDI"="CoMM"))
dat_rej$Method=droplevels(dat_rej$Method)
rho = 0.5; n2 = 8000;
t1e_rej = dat_rej[dat_rej$RhoX==rho&dat_rej$n2==n2,]

t1e_rej$h2z = factor(t1e_rej$h2z)
t1e_rej$h2y = factor(t1e_rej$h2y)
t1e_rej$Sparsity = factor(t1e_rej$Sparsity)
t1e_rej$n2 = factor(t1e_rej$n2)
t1e_rej$Method <- ordered(t1e_rej$Method, levels = c("CoMM","2-stage:Ridge","2-stage:Enet","SKAT"))
t1e_rej$Power = as.numeric(as.character((t1e_rej$Power)))

t1e_rej$h2y2 <- factor(t1e_rej$h2y, labels = c("h[C]^2==0.01", "h[C]^2==0.03", "h[C]^2==0.05", "h[C]^2==0.07", "h[C]^2==0.09"))
t1e_rej$h2z2 <- factor(t1e_rej$h2z, labels = c("h[T]^2==0", "h[T]^2==0.001", "h[T]^2==0.002", "h[T]^2==0.003"))

library(ggplot2)
ggplot(t1e_rej, aes(x = Sparsity, y = Power,fill = Method))+
geom_bar(stat="identity", position=position_dodge())+
geom_errorbar(aes(ymin=Power-sd_rej, ymax=Power+sd_rej), width=.2,
                 position=position_dodge(.9)) +
facet_grid(h2z2~h2y2,labeller = label_parsed,scales = "free_y") + 
geom_hline(yintercept=0.05,colour="orange",linetype="dashed")


