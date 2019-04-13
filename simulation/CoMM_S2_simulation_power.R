rm(list = ls())
library(SCLMM)
library(CoMM)
library(glmnet)
library(survey)
library(mvtnorm)
ls("package:CoMM")

setwd("/home/users/nus/gmsyan/SCLMM/File")
source("/home/users/nus/gmsyan/SCLMM/File/functionset.R")
L = 1; M = 100; #rho = 0.5;
n1 = 400; #n2 = 5000;
maf = as.numeric(unlist(read.table(paste("maf",M,".txt",sep=""))))
q1 = 1; q2 = 1; #h2y = 0.01;
nrep = 500;
rhoall = c(-0.8,-0.5,-0.2,0.2,0.5,0.8);
h2yall = c(0.01,0.03, 0.05, 0.07,0.09)
beta_propall = c(0.1,0.2,0.3,0.4,0.5,1)
h2all = c(0.001, 0.002, 0.003)
n2all = c(5000,8000)
n3all = c(400)
lam = 0.8
setwd("/home/users/nus/gmsyan/SCLMM/simulation/Result")

for (i in 4:6){
  for (j in 1:5){
    for (k in 1:6){
      for (l in 1){
        for (nn in 1){
          for(i_ in 1){
            
            rho = rhoall[i];
            h2y = h2yall[j];
            beta_prop = beta_propall[k];
            h2 = h2all[l];
            n2 = n2all[nn]
            n3 = n3all[i_]
            
            pval = matrix(nrow=nrep,ncol=12)
            alpha_hat = matrix(nrow=nrep,ncol=5)
            
            
            for (irep in 1:nrep){
              print(irep)
              X = genRawGeno(maf, L, M, rho, n1 + n2)
              X3 = genRawGeno(maf, L, M, rho, n3)
              
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
              sz2 <- var(y2*alpha)*((1 - h2)/h2)
              z <- alpha0 + y2*alpha + rnorm(n2,0,sqrt(sz2))
              
              y = y1;
              #normalize X1 and X2
              mean.x1 = apply(X1,2,mean);
              x1m = sweep(X1,2,mean.x1);
              std.x1 = apply(x1m,2,sd)
              x1p = sweep(x1m,2,std.x1,"/");
              #x1p = x1p/sqrt(dim(x1p)[2])
              
              mean.x2 = apply(X2,2,mean);
              x2m = sweep(X2,2,mean.x2);
              std.x2 = apply(x2m,2,sd)
              x2p = sweep(x2m,2,std.x2,"/");
              #x2p = x2p/sqrt(dim(x2p)[2])
              
              mean.x3 = apply(X3,2,mean);
              x3m = sweep(X3,2,mean.x3);
              std.x3 = apply(x3m,2,sd)
              x3p = sweep(x3m,2,std.x3,"/");
              #x3p = x3p/sqrt(dim(x3p)[2])
              
              w2 = matrix(rep(1,n2),ncol=1);
              w1 = matrix(rep(1,n1),ncol=1);
              
              ## CoMM
              fmHa = CoMM_covar_pxem(y, z, x1p, x2p, w1, w2, constr = 0, pxem_indicator = 0);
              fmH0 = CoMM_covar_pxem(y, z, x1p, x2p, w1, w2, constr = 1, pxem_indicator = 0);
              
              loglikHa = fmHa$loglik
              loglikH0 = fmH0$loglik
              
              tstat = 2 * (loglikHa - loglikH0);
              pval[irep,1] = pchisq(tstat,1,lower.tail=F)
              alpha_hat[irep,1] = fmHa$alpha
              
              pval[irep, 7] = fmHa$sigma2beta
              pval[irep, 8] = fmHa$sigma2y
              
              fit = rl_skat(x2p, z)
              pval[irep, 9] = fit$scores
              pval[irep, 10] = fit$pval
              
              hatmu1 = matrix(0, M, 1)
              hats1 = matrix(0, M, 1)
              
              for (m in 1:M){
                fm = lm(y~1+x1p[,m]);
                hatmu1[m] = summary(fm)$coefficients[2,1]
                hats1[m] = summary(fm)$coefficients[2,2];
              }
              
              opts = list(max_iter = 10000, dispF = 1, display_gap = 10, epsStopLogLik = 1e-5, fix_alphag = 0)
              opts1 = list(max_iter = 10000, dispF = 1, display_gap = 10, epsStopLogLik = 1e-5, fix_alphag = 1)
              
              w = matrix(0, n2, 1)
              objHa = SCLMM1(hatmu1, hats1, x2p, z, w, opts)
              objH0 = SCLMM1(hatmu1, hats1, x2p, z, w, opts1)
              pval[irep, 11] = objHa$Lq[length(objHa$Lq)] - objH0$Lq[length(objH0$Lq)]
              
              ## PrediXcan (ridge)
              cvfit1 = cv.glmnet(x1p, y, alpha = 0)
              beta_hat = coef(cvfit1, s = "lambda.min")
              yhat1 = predict(cvfit1,x2p,s="lambda.min")
              if (var(yhat1) != 0){
                fm_twas1 = summary(lm(z~cbind(yhat1,w2[,-1])));
                alpha_hat[irep,2] = fm_twas1$coefficients[2,1]
                pval[irep,2] = fm_twas1$coefficients[2,4]
              }
              if (var(yhat1) == 0){
                alpha_hat[irep,2] = NA;
                pval[irep,2] = NA;
              }
              
              ## PrediXcan (enet)
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
              
              sumx3p = apply(x3p*x3p, 2, sum)
              R = matrix(0, M, M);
              for (i1 in 1:M){
                for (j1 in 1:M){
                  R[i1,j1] = t(x3p[,i1])%*%x3p[,j1]/sqrt(sumx3p[i1]*sumx3p[j1])                        
                }
              }
              R = R*lam + (1 - lam)*diag(M)
              
              
              ## generate summary data
              hatmu = matrix(0, M, 1)
              hats = matrix(0, M, 1)
              
              for (m in 1:M){
                fm = lm(z~1+x2p[,m]);
                hatmu[m] = summary(fm)$coefficients[2,1]
                hats[m] = summary(fm)$coefficients[2,2];
              }
              
              ## CoMM-S^2
              opts = list(max_iter = 10000, dispF = 1, display_gap = 10, epsStopLogLik = 1e-5, fix_alphag = 0)
              opts1 = list(max_iter = 10000, dispF = 1, display_gap = 10, epsStopLogLik = 1e-5, fix_alphag = 1)
              
              px = 1
              objHa = SCLMM_S2(x1p, y, w1, hatmu, hats, R, opts, px)
              objH0 = SCLMM_S2(x1p, y, w1, hatmu, hats, R, opts1, px)
              
              stat = 2*(objHa$LRLB - objH0$LRLB)
              alpha_hat[irep,4] = objHa$alphag
              pval[irep, 4] = pchisq(stat, 1, lower.tail = F)
              
              
              ## S-PrediXcan
              varx3p = apply(x3p,2,var)
              if (var(x3p%*%beta_hat[-1]) != 0){
                Zg = sum(beta_hat[-1]*sqrt(varx3p)*hatmu/hats)/sqrt(var(x3p%*%beta_hat[-1]))
                pval[irep, 5] = pchisq(Zg^2, 1, lower.tail = F)
              }
              if (var(x3p%*%beta_hat[-1]) == 0){
                pval[irep,5] = NA;
              }
              
              varx3p = apply(x3p,2,var)
              if (var(x3p%*%beta_hat2[-1]) != 0){
                Zg = sum(beta_hat2[-1]*sqrt(varx3p)*hatmu/hats)/sqrt(var(x3p%*%beta_hat2[-1]))
                pval[irep, 6] = pchisq(Zg^2, 1, lower.tail = F)
              }
              if (var(x3p%*%beta_hat2[-1]) == 0){
                pval[irep,6] = NA;
              }
              
              
            }
            
            outfile1 <- paste("sim_power_SCLMM_pval_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                              q2,"_",h2y,"_",beta_prop,"_",h2,"_",lam, ".txt",sep="")
            outfile2 <- paste("sim_power_SCLMM_alpha_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                              q2,"_",h2y,"_",beta_prop,"_",h2,"_",lam, ".txt",sep="")
            write.table(pval,outfile1,sep="\t",quote=F, row.names=F,col.names=F)
            write.table(alpha_hat,outfile2,sep="\t",quote=F, row.names=F,col.names=F)
            
          }
        }
      }
      
    }
    
  }
  
}
