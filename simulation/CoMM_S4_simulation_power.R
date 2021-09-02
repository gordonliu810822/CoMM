rm(list = ls())
library(CoMM)
library(glmnet)
library(survey)
library(mvtnorm)
library(mammot)
library(Rcpp)
ls("package:CoMM")


setwd("/home/projects/11000368/Work/YangYi/SCLMM/File")
source("/home/projects/11000368/Work/YangYi/SCLMM/File/functionset.R")
L = 1; M = 100; #rho = 0.5;
n1 = 5000; #n2 = 5000;
maf = as.numeric(unlist(read.table(paste("maf",M,".txt",sep=""))))
q1 = 1; q2 = 1; #h2y = 0.01;
nrep = 2000;
rhoall = c(-0.8,-0.5,-0.2,0.2,0.5,0.8);
h2yall = c(0.01,0.03, 0.05, 0.07,0.09)
beta_propall = c(0.1,0.2,0.3,0.4,0.5,1)
h2all = c(0.001, 0.002, 0.003)
n2all = c(5000,8000)
n3all = c(400)
lam = 0.8
setwd("/home/projects/11000368/Work/YangYi/SCLMM/simulation/Result")

for (i in 4){
  for (j in 1:5){
    for (k in c(2,6)){
      for (l in 1){
        for (nn in 1){
          for(i_ in 1){
            i = 1; j = 1; k = 1; l = 1; nn = 1; i_ = 1;
            
            rho = rhoall[i];
            h2y = h2yall[j];
            beta_prop = beta_propall[k];
            h2 = h2all[l];
            n2 = n2all[nn]
            n3 = n3all[i_]
            
            pval = matrix(nrow=nrep,ncol=2)
            alpha_hat = matrix(nrow=nrep,ncol=2)
            
            
            for (irep in 1:nrep){
              print(irep)
              X = genRawGeno(maf, L, M, rho, n1 + n2)
              
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
              alpha0 <- 0
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
              
              X3 = X1[1:400,]
              mean.x3 = apply(X3,2,mean);
              x3m = sweep(X3,2,mean.x3);
              std.x3 = apply(x3m,2,sd)
              x3p = sweep(x3m,2,std.x3,"/");
              #x3p = x3p/sqrt(dim(x3p)[2])
              
              X4 = X2[1:400,]
              mean.x4 = apply(X4,2,mean);
              x4m = sweep(X4,2,mean.x4);
              std.x4 = apply(x4m,2,sd)
              x4p = sweep(x4m,2,std.x4,"/");
              
              w2 = matrix(rep(1,n2),ncol=1);
              w1 = matrix(rep(1,n1),ncol=1);
              
              ## R1
              sumx3p = apply(x3p*x3p, 2, sum)
              R = matrix(0, M, M);
              for (i1 in 1:M){
                for (j1 in 1:M){
                  R[i1,j1] = t(x3p[,i1])%*%x3p[,j1]/sqrt(sumx3p[i1]*sumx3p[j1])                        
                }
              }
              R = R*lam + (1 - lam)*diag(M)
              
              ## R2
              sumx4p = apply(x4p*x4p, 2, sum)
              R2 = matrix(0, M, M);
              for (i1 in 1:M){
                for (j1 in 1:M){
                  R2[i1,j1] = t(x4p[,i1])%*%x4p[,j1]/sqrt(sumx4p[i1]*sumx4p[j1])                        
                }
              }
              R2 = R2*lam + (1 - lam)*diag(M)
              
              
              ## generate summary data
              ## summary statistics for eQTL data
              hatmu = matrix(0, M, 1)
              hats = matrix(0, M, 1)
              
              for (m in 1:M){
                fm = lm(y~1+x1p[,m]);
                hatmu[m] = summary(fm)$coefficients[2,1]
                hats[m] = summary(fm)$coefficients[2,2];
              }
              
              ## summary statistics for GWAS data
              hatmu2 = matrix(0, M, 1)
              hats2 = matrix(0, M, 1)
              
              for (m in 1:M){
                fm = lm(z~1+x2p[,m]);
                hatmu2[m] = summary(fm)$coefficients[2,1]
                hats2[m] = summary(fm)$coefficients[2,2];
              }
              
             
              ## CoMM-S^4
              opts = list(max_iter = 10000, dispF = 1, display_gap = 10, epsStopLogLik = 1e-5, fix_alphag = 0)
              opts1 = list(max_iter = 10000, dispF = 1, display_gap = 10, epsStopLogLik = 1e-5, fix_alphag = 1)
              
              score = hatmu/hats
              score2 = hatmu2/hats2
              objHa = CoMM_S4(score, score2, R, R2, opts)
              objH0 = CoMM_S4(score, score2, R, R2, opts1)
              
              stat2 = 2*(objHa$LRLB - objH0$LRLB)
              alpha_hat[irep,2] = objHa$alphag
              pval[irep, 2] = pchisq(stat2, 1, lower.tail = F)
              
            }
            
            outfile1 <- paste("sim_scatter_power_CoMM_2S_pval_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                              q2,"_",h2y,"_",beta_prop,"_",h2,"_",lam, ".txt",sep="")
            outfile2 <- paste("sim_scatter_power_CoMM_2S_alpha_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                              q2,"_",h2y,"_",beta_prop,"_",h2,"_",lam, ".txt",sep="")
            write.table(pval,outfile1,sep="\t",quote=F, row.names=F,col.names=F)
            write.table(alpha_hat,outfile2,sep="\t",quote=F, row.names=F,col.names=F)
            
          }
        }
      }
      
    }
    
  }
  
}
