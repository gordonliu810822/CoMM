# kernel base rl_skat, X is fixed covariates, Z is random covariates, K = Z %*% W %*% t(Z)
require(survey)

#' @title
#' CoMM
#' @description
#' fit SKAT using RL-SKAT: An Exact and Efficient Score Test for Heritability and Set Tests.
#'
#' @param y  a vector for the phenotype of GWAS.
#' @param Z  a standardized genotype matrix for GWAS data.
#' @param Z  a matrix of coveriates.
#'
#' @return List of model parameters: test statistics (scores) and p-value (pval)
#'
#' @export
rl_skat <- function(Z,y,X=NULL,intercept=T,method=NULL){
  y <- as.matrix(y)
  N <- ncol(y)

  n <- nrow(Z)
  m <- ncol(Z)
  if(intercept){
    if(is.null(X)){
      X <- matrix(1,nrow = n,ncol = 1)
    } else {
      X <- cbind(rep(1,n),X)
    }
  }

  if(is.null(X)){
    p <- 0
  } else {
    p <- ncol(X)
    pinvX <- solve(t(X)%*%X) %*% t(X)
  }

  if(is.null(method)){
    if(n > m){
      method <- "low-rank"
    } else {
      method <- "kernel"
    }
  }

  
  
  if(method=="kernel"){
    K <- Z %*% t(Z)
    # calculate SKS
    if(is.null(X)){
      SKS <- K
    } else {
      SKS <- K - X %*% pinvX %*% K
      SKS <- SKS - SKS %*% X %*% pinvX
    }
    
    #Perform eigen decomposition
    phi <- eigen(SKS,only.values = T)$values
    
    phi[phi < 1e-5] <- 0
    phi <- sort(phi,decreasing = T)
    
    #calculate rank
    k <- sum(phi > 0)
    
    #calculate q
    if(is.null(X)){
      q <- n - k
    } else{
      B <- cbind(K,X)
      q <- n - sum(svd(B)$d > 1e-5)
    }
    
    #calculate scores
    numerator <- colSums(SKS %*% y * y)
    if(is.null(X)){
      denominator <- colSums(y^2)
    } else {
      denominator <- colSums((y-X%*%pinvX%*%y) * y)
    }
    scores <- numerator / denominator * (n-p)
  }

  if(method=="low-rank"){
    #calculate SZ
    if(is.null(X)){
      SZ <- Z
    } else {
      SZ <- Z - X %*% (pinvX %*% Z)
    }
    
    #Perform eigen decomposition
    phi <- svd(SZ)$d^2
    phi <- c(phi,rep(0,n-length(phi)))
    
    phi[phi < 1e-5] <- 0
    phi <- sort(phi,decreasing = T)
    
    #calculate rank
    k <- sum(phi > 0)
    
    #calculate q
    if(is.null(X)){
      q <- n - k
    } else {
      B <- cbind(Z,X)
      q <- n - sum(svd(B)$d > 1e-5)
    }
    
    #calculate scores
    if(is.null(X)) {
      y_proj <- y
    } else {
      y_proj <- y - X %*% (pinvX %*% y)
    }
    numerator <- colSums((t(Z) %*% y_proj)^2)
    denominator <- colSums(y_proj^2)
    scores <- numerator / denominator * (n-p)
  }
  

  #calculate p-values
  pval <- sapply(scores,
                 function(r,phi,n,p,q) pchisqsum(0,df=rep(1,k+q),a=c(phi[1:k]-r/(n-p),-rep(r/(n-p),q)),
                                                 method="integration",lower.tail = F),
                 phi=phi,n=n,p=p,q=q)

  ret <- list(scores=scores,pval=pval)
}
