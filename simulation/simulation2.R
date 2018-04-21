#set.seed(1000);
#setwd("D:\\Data\\GoogleDrive\\Work\\Multiple_Trait\\TWAS\\AUDI\\Simulation")
#for (nsnp in c(100,200,400)){
#    maf = runif(nsnp,0.05,0.5)
#    savefile3=paste("maf",nsnp,".txt",sep="")
#    write.table(maf,savefile3,sep="\t",quote=F, row.names=F,col.names=F)
#}


rm(list = ls())
library(AUDI)
library(glmnet)
library(survey)
ls("package:AUDI")
#setwd("G:\\Data\\GoogleDrive\\Work\\Multiple_Trait\\TWAS\\AUDI\\Simulation")
#setwd("/scratch/users/nus/gmsjl/AUDI/sim")
setwd("/home/projects/11000369/Work/AUDI/sim")
source("functionset.R")
L = 1; M = 100; #rho = 0.5;
n1 = 350; #n2 = 5000;
maf = as.numeric(unlist(read.table(paste("maf",M,".txt",sep=""))))
q1 = 1; q2 = 1; #h2y = 0.01;
#beta_prop = 0.2;
#alpha_true= 0.2;
nrep = 500;
rhoall = c(0.2,0.5,0.8);
h2yall = seq(0.005,0.1,by=0.005)
beta_propall = c(0.2,0.4,1)
h2all = c(0.001, 0.002, 0.003);
n2all = c(4000,8000)
#setwd("/scratch/users/nus/gmsjl/AUDI/sim/results")
setwd("/home/projects/11000369/Work/AUDI/sim/sim2")
#sim_AUDI_pval_350_10000_1_100_0.2_1_1_0.25_0.01_0.001.txt
# i =1; j =1; k=1; l = 3; nn = 2;
for ( i in 1:3){
    for (j in 1){
        for (k in 1){
            for ( l in 1){
               for (nn in 1:2){

                rho = rhoall[i];
                h2y = h2yall[j];
                beta_prop = beta_propall[k];
                h2 = h2all[l];
                n2 = n2all[nn]

                pval = matrix(nrow=nrep,ncol=5)
                alpha_hat = matrix(nrow=nrep,ncol=4)


                for (irep in 1:nrep){
                  print(irep)
                  X = genRawGeno(maf, L, M, rho, n1 + n2)

                  #X1 = X[1:n1,];
                  #X2 = X[(n1+1):(n1+n2),]

                  #dat = genPheno(X,n1,n2, q1, q2, h2y, beta_prop, h2)
                  # var(dat$x1p%*%dat$beta)/var(dat$y)
                  # var(dat$y02)/var(dat$z)
                  #z = dat$z;
                  #y = dat$y;
                  #w1 = dat$w1;
                  #w2 = dat$w2;
                  #x1p = dat$x1p;
                  #x2p = dat$x2p;
                  #SIGMA = matrix(nrow=M,ncol=M)
                  #for (ii in 1:M){
                  #    for (jj in 1:M){
                  #        SIGMA[ii,jj] = rho^(abs(ii-jj));
                  #    }
                  #}
                  #X = rmvnorm(n1 + n2, mean=rep(0,M), sigma=SIGMA, method="chol")
                  b = numeric(M);
                  m = M * beta_prop;
                  b[sample(M,m)] = rnorm(m);

                  b0 = 6
                  y0 <- X%*%b + b0
                  y  <- y0 + (as.vector(var(y0)*(1-h2y)/h2y))^0.5*rnorm(n1+n2)

                  y1 <- y[1:n1]
                  X1 <- X[1:n1,]
                  y2 <- y0[(n1+1):(n1+n2)]
                  X2 <- X[(n1+1):(n1+n2),]
                  alpha0 <- 3
                  alpha <- 0.3
                  sz2 <- var(y2*alpha) * ((1-h2)/h2)
                  z <- alpha0 + y2*alpha + rnorm(n2,0,sqrt(sz2))

                  #x1p=X1; x2p = X2;
                  y = y1;
                  #normalize X1 and X2
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

                  ## AUDI
                  fm0 = lmm_pxem(y, w1,x1p, 100)
                  sigma2beta =fm0$sigma2beta;
                  sigma2y =fm0$sigma2y;
                  beta0 = fm0$beta0;

                  fmHa = AUDI_covar_pxem(y, z, x1p, x2p, w1, w2, sigma2beta, sigma2y, beta0, 0, 1e-5, 1000);
                  fmH0 = AUDI_covar_pxem(y, z, x1p, x2p, w1, w2, sigma2beta, sigma2y, beta0, 1, 1e-5, 1000);

                  loglikHa = max(fmHa$loglik,na.rm=T)
                  loglikH0 = max(fmH0$loglik,na.rm=T)
                  #tstat[irep,1] = 2 * (loglikHa - loglikH0);
                  tstat = 2 * (loglikHa - loglikH0);
                  pval[irep,1] = pchisq(tstat,1,lower.tail=F)
                  alpha_hat[irep,1] = fmHa$alpha

                  ## 2stage AUDI
                  fm2Ha = covar_pxem_2ndstage(y,z,x1p, x2p, w1, w2, fm0$sigma2beta, fm0$sigma2y, fm0$beta0,0,1e-6,20);
                  fm2H0 = covar_pxem_2ndstage(y,z,x1p, x2p, w1, w2, fm0$sigma2beta, fm0$sigma2y, fm0$beta0,1,1e-6,20);
                  tstat2 = 2*(max(fm2Ha$loglik,na.rm=T)-max(fm2H0$loglik,na.rm=T))
                  pval[irep,5] = pchisq(tstat2,1,lower.tail=F)
                  alpha_hat[irep,4] = fm2Ha$alpha

                  ## TWAS (ridge)
                  cvfit1 = cv.glmnet(x1p, y, alpha = 0)
                  #beta_hat1 = coef(cvfit1, s = "lambda.min")
                  yhat1 = predict(cvfit1,x2p,s="lambda.min")
                  #yhat1 = x2p%*%beta_hat1[-1] + beta_hat1[1];
                  if (var(yhat1) != 0){
                     fm_twas1 = summary(lm(z~cbind(yhat1,w2[,-1])));
                     alpha_hat[irep,2] = fm_twas1$coefficients[2,1]
                     #tstat[irep,2] = fm_twas$coefficients[2,3]^2
                     pval[irep,2] = fm_twas1$coefficients[2,4]
                  }
                  if (var(yhat1) == 0){
                     alpha_hat[irep,2] = NA;
                     pval[irep,2] = NA;
                  }

                  ## TWAS (enet)
                  cvfit2 = cv.glmnet(x1p, y, alpha = 0.5)
                  beta_hat2 = coef(cvfit2, s = "lambda.min")

                  yhat2 = x2p%*%beta_hat2[-1] + beta_hat2[1];
                  if (var(yhat2) != 0){
                     fm_twas2 = summary(lm(z~cbind(yhat2,w2[,-1])));
                     alpha_hat[irep,3] = fm_twas2$coefficients[2,1]
                     pval[irep,3] = fm_twas2$coefficients[2,4]
                  }
                  if (var(yhat2) == 0){
                     alpha_hat[irep,3] = NA;
                     pval[irep,3] = NA;
                  }
                  ## SKAT
                  fm_skat <- rl_skat(X2,z)
                  pval[irep,4] = fm_skat$pval
                }

                outfile1 <- paste("sim2_AUDI_pval_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")
                outfile2 <- paste("sim2_AUDI_alpha_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")
                write.table(pval,outfile1,sep="\t",quote=F, row.names=F,col.names=F)
                write.table(alpha_hat,outfile2,sep="\t",quote=F, row.names=F,col.names=F)

               }
            }

        }

    }

}
#apply(alpha_hat,2,mean)
#boxplot(alpha_hat)
#apply(tstat >qchisq(0.95,1),2,sum)
