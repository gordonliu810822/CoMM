rm(list = ls())
library(CoMM)
library(glmnet)
library(survey)
library(mvtnorm)
ls("package:CoMM")

L = 1; M = 100; #rho = 0.5;
n1 = 350; #n2 = 5000;
maf = as.numeric(unlist(read.table(paste("maf",M,".txt",sep=""))))
q1 = 1; q2 = 1; #h2y = 0.01;
nrep = 1000;
rhoall = c(-0.8,-0.5,-0.2,0.2,0.5,0.8);
h2yall = c(0.01,0.03, 0.05, 0.07,0.09)
beta_propall = c(0.1,0.2,0.3,0.4,0.5,1)
h2all = c(0)
n2all = c(4000,8000)

for ( i in 1:6){
    for (j in 1){
        for (k in 1){
            for ( l in 1){
               for (nn in 1){

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
                  alpha <- 0
                  sz2 <- var(y2)
                  z <- alpha0 + y2*alpha + rnorm(n2,0,sqrt(sz2))

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

                  ## CoMM
                  fmHa = CoMM_covar_pxem(y, z, x1p, x2p, w1, w2, constr = 0);
                  fmH0 = CoMM_covar_pxem(y, z, x1p, x2p, w1, w2, constr = 1);

                  loglikHa = fmHa$loglik
                  loglikH0 = fmH0$loglik

                  tstat = 2 * (loglikHa - loglikH0);
                  pval[irep,1] = pchisq(tstat,1,lower.tail=F)
                  alpha_hat[irep,1] = fmHa$alpha
                  
                  ## PrediXcan (ridge)
                  cvfit1 = cv.glmnet(x1p, y, alpha = 0)
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
                  ## SKAT
                  fm_skat <- rl_skat(X2,z)
                  pval[irep,4] = fm_skat$pval
                }

                outfile1 <- paste("sim_t1e_CoMM_pval_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")
                outfile2 <- paste("sim_t1e_CoMM_alpha_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
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
