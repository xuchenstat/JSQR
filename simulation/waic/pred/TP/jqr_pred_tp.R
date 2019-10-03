tau.grid <- c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99)
tau.len <- length(tau.grid)

n.test <- 50

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

corfun <- function(phi, distance, kernel, kappa, cte){
  if(kernel == 1){
    corEntry <- sapply(distance, sqexp, phi = phi)
  } else {
    corEntry <- sapply(distance, matern, phi = phi, kappa = kappa, cte = cte)
  }
  return(matrix(corEntry, nrow = 1))
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


nphi <- 10

library(kernlab)
library(Matrix)
library(qrjoint)
library(mvtnorm)


jqr.test.loss <- jqr.test.qt <- array(NA, dim = c(100, n.test, tau.len))


for(datanum in 601:700){
  
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
  
  z <- t(chol((alpha*C+(1-alpha)*diag(n))/rgamma(1, 3/2, 3/2)))%*%rnorm(n)
  tau <- pt(z, df = 3)
  x <- matrix(runif(n*p,-1,1), nrow = n)
  y <- rep(NA, n)
  for(i in 1:n){
    b0 <- integrate(b0.prime, 0.5, tau[i], rel.tol = 1e-10)$value + b0.half
    b <- integrate(b.prime, lower = 0.5, upper = tau[i], rel.tol = 1e-10)$value + b.half
    y[i] <- b0 + b*x[i]
  }
  
  location.test <- matrix(runif(2*n.test), nrow = n.test)
  x.test <- matrix(runif(n.test*p,-1,1), nrow = n.test)
  
  distance.test <- t(apply(location.test, 1, function(u) apply(location, 1, function(v) sqrt(sum((u-v)^2)))))
  
  sigma22 <- alpha*C+(1-alpha)*diag(n)
  sigma22.inv <- chol2inv(chol(sigma22))
  dmaha <- t(z)%*%sigma22.inv%*%z
  z.grid <- qt(tau.grid, df = 3 + n)
  
  tau.test <- qt.test <- matrix(NA, nrow = n.test, ncol = tau.len)
  
  for(i in 1:n.test){
    sigma11 <- 1
    sigma12 <- alpha*corfun(phi = phi, distance = distance.test[i,], kernel = 2, kappa = kappa, cte = cte)
    sigma12.22.inv <- sigma12%*%sigma22.inv
    cond.var <- as.numeric((3 + dmaha)*(sigma11 - tcrossprod(sigma12.22.inv, sigma12))/(3 + n))
    cond.mean <- c(sigma12.22.inv %*% z)
    tau.test[i,] <- pt(cond.mean + sqrt(cond.var) * z.grid, df = 3)
    for(k in 1:tau.len){
      b0 <- integrate(b0.prime, 0.5, tau.test[i,k], rel.tol = 1e-10)$value + b0.half
      b <- integrate(b.prime, lower = 0.5, upper = tau.test[i,k], rel.tol = 1e-10)$value + b.half
      qt.test[i,k] <- b0 + b*x.test[i]
    }
  }
  
  
  
  ######################## jqr
  file.name <- paste('jqrfit', datanum, '.RData', sep = '')
  load(file.name)
  
  coef.jqr <- coef(jqr.fit, plot = FALSE, nmc = 500)
  coef.jqr.mean <- apply(coef.jqr$beta.samp[as.integer(tau.grid/0.01)+1,,], c(1,2), mean)
  
  jqr.test.qt[datanum - 600,,] <- apply(coef.jqr.mean, 1, function(u) as.matrix(cbind(1,x.test))%*%u)
  jqr.test.loss[datanum - 600,,] <- abs(jqr.test.qt[datanum - 600,,] - qt.test)
  
  print(datanum)
}