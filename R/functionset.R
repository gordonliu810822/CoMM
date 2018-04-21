## simulation maf first
library(mvtnorm)
genRawGeno <- function(maf, L, M, rho, n){
   SIGMA = matrix(nrow=M,ncol=M)
   for (i in 1:M){
       for (j in 1:M){
           SIGMA[i,j] = rho^(abs(i-j));
       }
   }

   nsnp = L*M;
   X = NULL;
   for ( l in 1:L ){
       #set.seed(1000);
       index = (M*(l-1)+1): (M*l);
       AAprob = maf[index]^2.;
       Aaprob = 2*maf[index]*(1-maf[index]);
       quanti = matrix(c(1-Aaprob-AAprob, 1- AAprob),M,2);
       Xt = rmvnorm(n, mean=rep(0,M), sigma=SIGMA, method="chol")
       Xt2 = matrix(0,n,M);
       for (j in 1:M){
           cutoff = qnorm(quanti[j,]);
           Xt2[Xt[,j] < cutoff[1],j] = 0;
           Xt2[Xt[,j] >= cutoff[1] & Xt[,j] < cutoff[2],j] = 1;  ## attention
           Xt2[Xt[,j] >= cutoff[2],j] = 2;
       }
       X <- cbind(X,Xt2);
   }
   return(X)
}

genPheno <- function(X,n1,n2, q1, q2, h2y, beta_prop, h2){
   X1 = X[1:n1,];
   X2 = X[(n1+1):(n1+n2),]
   n1 = dim(X1)[1];
   n2 = dim(X2)[1];
   p1 = dim(X1)[2];
   p2 = dim(X2)[2];

   if (p1 != p2){
      err("Dimension in X1 and X2 are not matched!")
   }

   p = p1;

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

   if ( q1 > 1){
      w1 = cbind(rep(1,n1),matrix(rnorm(n1*(q1-1)),nrow=n1,ncol=q1-1));
   }
   else if (q1 == 1){
      w1 = matrix(rep(1,n1),ncol=1);
   }

   if ( q2 > 1){
      w2 = cbind(rep(1,n2),matrix(rnorm(n2*(q2-1)),nrow=n2,ncol=q2-1));
   }
   else if (q2 == 1){
      w2 = matrix(rep(1,n2),ncol=1);
   }

   alpha0 = rnorm(q2);
   beta0 = 2*rnorm(q1);

   m = p * beta_prop;
   beta = numeric(p);
   beta[sample(p,m)] = rnorm(m);
   #ey = matrix(rnorm(n1+n2),ncol=1);
   #y0 = x1p%*%beta + w1%*%beta0;
   #y = y0 + (ey[1:n1])%*%unlist(sqrt(var(x1p%*%beta))*sqrt((1-h2y)/h2y));

   #ez = rnorm(n2);
   #y02 = x2p%*%beta + w2%*%alpha0 +(ey[(n1+1):(n1+n2)])%*%unlist(sqrt(var(x1p%*%beta))*sqrt((1-h2y)/h2y));
   y0 <- X%*%b + b0
   y  <- y0 + (as.vector(var(y0)*(1-h2y)/h2y))^0.5*rnorm(n1+n2)
   y1 <- y[1:n1]
   y2 <- y0[(n1+1):(n1+n2)]

   alpha0 <- 3
   alpha <- 0.3



   if (h2 !=0 ){
      sz2 <- var(y2*alpha) * ((1-h2)/h2)
      z <- alpha0 + y2*alpha + rnorm(n2,0,sqrt(sz2))
   } 
   else if (h2 == 0){
     z = rnorm(n2);
   }
   y = y1;

   return(list(beta=beta,z=z,y=y,x1p=x1p,x2p=x2p,w1=w1,w2=w2))
}
