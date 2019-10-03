library(ald)
library(plyr)
library(quantreg)
library(splines)
library(kernlab)
library(Matrix)


dyn.load('/jsqr/jsqrtp/jsqrtp.so')
source('/jsqr/jsqrtp/utility.R')
source('/jsqr/jsqrtp/jsqrtp.R')


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


b0.half <- 0
b.half <- 0

b0.prime <- function(tau){3/dt(qt(tau, df = 3), df = 3)}

v <- function(tau){
  vtau <- 3*(tau - 0.5)
  return(vtau/sqrt(1+vtau^2))
}


b.prime <- function(tau){
  return(b0.prime(tau)*v(tau))
}



for (datanum in 501:600){
  
  set.seed(datanum)    
  
  n <- 500; p <- 1; location <- matrix(runif(2*n), nrow = n)
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
  
  z <- t(chol(alpha*C+(1-alpha)*diag(n)))%*%rnorm(n)
  tau <- pnorm(z)
  x <- matrix(runif(n*p,-1,1), nrow = n)
  y <- rep(NA, n)
  for(i in 1:n){
    b0 <- integrate(b0.prime, 0.5, tau[i], rel.tol = 1e-10)$value + b0.half
    b <- integrate(b.prime, lower = 0.5, upper = tau[i], rel.tol = 1e-10)$value + b.half
    y[i] <- b0 + b*x[i]
  }  
  
  jsqr.fit <- qrjointtp(x, y, par = 'prior', distance = distance, nphi = 10, kernel = "matern", acpt.hptarget = 0.15, acpt.target = 0.15,  kappa = 2, nsamp = 2e3, thin = 10, fbase = 't')

  save(jsqr.fit, file = paste('jsqrtpfit', datanum, '.RData', sep = ''))
}