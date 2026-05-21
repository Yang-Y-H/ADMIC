# DATA is an abundance matrix where the rows correspond to taxa and the columns correspond to samples. 
# rho is (Non-negative) regularization parameter for lasso. The default is 0.1.
library(glasso)
bonobo.precision <- function(DATA, rho = 0.1){
  n.sample <- ncol(DATA)
  g <- nrow(DATA)
  DATA[which(DATA == 0)] <- 1
  x <- log(DATA)
  x.mean <- rowMeans(x)
  X <- x-x.mean
  leaveoneCOV <- function(i){
    data <- t(X[,-i])
    return(cov(data))
  }
  S <- lapply(1:n.sample, leaveoneCOV)
  diagS <- matrix(0, nrow = g, ncol = n.sample)
  for(i in 1:n.sample){
    diagS[,i] <- diag(S[[i]])
  }
  eta <- rep(0, g)
  for(j in 1:g){
    eta[j] <- var(diagS[j,])*(n.sample-1)/n.sample
  }
  eta.sum <- sum(eta)
  sk <- rep(0, n.sample)
  for(i in 1:n.sample){
    sk[i] <- 2*sum(diagS[,i]^2)
  }
  vi_g_3 <- sk/eta.sum
  vi_g <- vi_g_3+3
  delta <- 1/vi_g
  
  sample.specific.cov <- function(i){
    Vi <- delta[i]*X[,i]%*%t(X[,i])+(1-delta[i])*S[[i]]
    diag(Vi)[which(diag(Vi)==0)] <- 1
    return(Vi)
  }
  V <- lapply(1:n.sample, sample.specific.cov)
  
  sample.specific.precision <- function(i){
    precision <- glasso(V[[i]], rho=rho)[["wi"]]
    return(precision)
  }
  Precision <- lapply(1:n.sample, sample.specific.precision)
  
  return(Precision)
}
