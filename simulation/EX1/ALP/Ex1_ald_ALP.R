library(ald)
library(plyr)
library(splines)
library(kernlab)
library(Matrix)
library(bayesQR)


covmat <- function(phi, n, distance, kernel, kappa, cte){
  cov.mat <- diag(n)
  if(kernel == 1){
    corEntry <- sapply(distance, sqexp, phi = phi)
  } else {
    corEntry <- sapply(distance, matern, phi = phi, kappa = kappa, cte = cte)
  }
  s <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n) {
      cov.mat[i,j] <- corEntry[s]
      cov.mat[j,i] <- cov.mat[i,j]
      s <- s + 1
    }
  }
  return(cov.mat)
}


matern <- function(x, phi, kappa, cte){
  uphi <- sqrt(2*kappa)*x/phi;
  if(uphi == 0.0) {
    return(1)
  } else {
    if(kappa == 0.5) {
      ans <- exp(-uphi)
    } else {
      ans <- cte * uphi^kappa * besselK(uphi, kappa, FALSE)
    }
  }
  return(ans)
}

tau.grid <- c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99)

for (datanum in 1:100){
  
  set.seed(datanum)    
  
  p.ald <- 0.4
  tau.ald <- 1
  
  n <- 200; p <- 1; location <- matrix(runif(2*n), nrow = n)
  
  distance = as.vector(dist(location, method = 'euclidean'))
  phi <- 0.3; kappa <- 2
  cte <- 2^(1-kappa)/gamma(kappa)
  corEntry <- sapply(distance, matern, phi = phi, kappa = kappa, cte = cte)
  C <- diag(n)
  s <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n) {
      C[i,j] <- corEntry[s]
      C[j,i] <- C[i,j]
      s <- s + 1
    }
  }
  
  alpha <- 0.7
  
  z.ald <- t(chol(alpha*C+(1-alpha)*diag(n)))%*%rnorm(n)
  xi.ald <- rexp(n, tau.ald)
  tau <- pALD(sqrt(2*xi.ald/tau.ald/p.ald/(1-p.ald))*z.ald+(1-2*p.ald)/p.ald/(1-p.ald)*xi.ald, mu = 0, sigma = 1/tau.ald)
  
  
  x <- matrix(runif(n*p,-1,1), nrow = n)
  y <- -3*(tau-0.5)*log(tau*(1-tau))-4*(tau-0.5)^2*log(tau*(1-tau))*x
  
  df <- data.frame(cbind(y,x))
  
  ald.fit <- bayesQR(formula = X1 ~., quantile = tau.grid, data = df, normal.approx = FALSE, ndraw = 2e4, keep = 20)
  save(ald.fit, file = paste('aldfit', datanum, '.RData', sep = ''))
}