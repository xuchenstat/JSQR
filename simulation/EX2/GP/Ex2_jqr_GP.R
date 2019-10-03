library(ald)
library(plyr)
library(quantreg)
library(splines)
library(kernlab)
library(Matrix)
library(qrjoint)

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
b.half <- c(0.96,-0.38,0.05,-0.22,-0.80,-0.80,-5.97)

a <- matrix(c(0,-3,0,0,0,-2,-3,0,2,-2,2,2,0,4,-4,5,1,0,-1,0,0), nrow = 3)

b0.prime <- function(tau){1/dnorm(qnorm(tau))}

v <- function(tau, j){
  vtau <- sapply(tau, function(u) t(a)%*%dnorm(u, c(0,0.5,1), 1/3))
  return(vtau[j,]/sqrt(1+apply(vtau^2,2,sum)))
}


b.prime <- function(tau, j){
  return(b0.prime(tau)*v(tau, j))
}


for (datanum in 301:400){
  
  set.seed(datanum)    
  
  n <- 500; p <- 7; location <- matrix(runif(2*n), nrow = n)
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
  
  z1 <- matrix(rnorm(n*p),n,p)
  zn <- sweep(z1,1,sqrt(apply(z1^2,1,sum)),'/')
  u <- rbeta(n,p,1)
  x <- sweep(zn,1,u,'*')
  
  y <- rep(NA, n)
  for(i in 1:n){
    b0 <- integrate(b0.prime, 0.5, tau[i], rel.tol = 1e-10)$value + b0.half
    b <- sapply(1:p, function(u) integrate(b.prime, lower = 0.5, upper = tau[i], j = u, rel.tol = 1e-10)$value) + b.half
    y[i] <- b0 + b%*%x[i,]
  }
  
  df <- data.frame(cbind(y,x))
  jqr.fit <- qrjoint(formula = y~., data = df, nsamp = 2e3, fbase = 'logistic')
  
  save(jqr.fit, file = paste('jqrfit', datanum, '.RData', sep = ''))
}