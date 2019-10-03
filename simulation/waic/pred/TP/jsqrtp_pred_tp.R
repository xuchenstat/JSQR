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
library(mvtnorm)
library(doParallel)
library(foreach)

source('/jsqr/jsqrtp/utility.R')

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
  
  
  ######################## jsqrtp
  file.name <- paste('jsqrtpfit', datanum, '.RData', sep = '')
  load(file.name)
  
  nsamp <- jsqr.fit$dim[11]; burn.perc <- 0.5; nmc <- 5e2; m <- jsqr.fit$dim[5]
  usamp <- matrix(jsqr.fit$usamp,ncol=n,byrow=T)
  parmgsamp <- matrix(jsqr.fit$parmgsamp, nrow = nsamp, byrow = TRUE)
  coefgp <- coef.qrjointgp(jsqr.fit, plot=FALSE, nmc = nmc, burn.perc = burn.perc, reduce = FALSE)
  ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 1)[-1]))
  
  phi.lb <- uniroot(function(u) sapply(u, matern, x = max(distance)/4, kappa = kappa, cte = cte) - 0.05, interval = c(1e-4, 10*max(distance)), tol = 1e-10)$root
  phi.ub <- 3*phi.lb
  phi.grid <- seq(phi.lb, phi.ub, length.out = nphi)
  
  tau.g <- jsqr.fit$tau.g
  len.reg.ix <- length(tau.g)
  tau.post <- array(NA, dim = c(nmc, n.test, len.reg.ix))
  tau.train.post <- array(NA, dim = c(nmc, n, len.reg.ix))
  
  cov.grid <- array(NA, dim = c(nphi, n, n))
  for(i in 1:nphi){
    cov.grid[i,,] <- covmat(phi = phi.grid[i], distance = distance, n = n, kernel = 2, kappa = kappa, cte = cte)
  }
  
  for(i in 1:nmc){
    psi <- plogis(jsqr.fit$lpsisamp[ss[i]])*18+2
    z.g <- qt(tau.g, df = psi)
    zsamp <- qt(usamp[ss[i],], df = psi)
    
    alpha <- exp(jsqr.fit$lkappasamp[ss[i]])/exp(parmgsamp[ss[i],(m+1)*(p+1)+1])
    sigma22 <- alpha*cov.grid[jsqr.fit$phisamp[ss[i]]+1,,]+(1-alpha)*diag(n)
    sigma22.inv <- chol2inv(chol(sigma22))
    dmaha <- t(zsamp)%*%sigma22.inv%*%zsamp
    
    for(j in 1:n.test){
      sigma11 <- 1
      sigma12 <- alpha*corfun(phi = phi.grid[jsqr.fit$phisamp[ss[i]]+1], distance = distance.test[j,], kernel = 2, kappa = kappa, cte = cte)
      sigma12.22.inv <- sigma12%*%sigma22.inv
      cond.var <- as.numeric((psi + dmaha)*(sigma11 - tcrossprod(sigma12.22.inv, sigma12))/(psi + n))
      cond.mean <- c(sigma12.22.inv %*% zsamp)
      tau.post[i,j,] <- pt((z.g - cond.mean)/sqrt(cond.var), df = psi + n)
    }
  }
  rm(cov.grid)
  gc()
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  betahat <- function(k, n.loc, n.col, len.tau.comp){
    beta.hat <- array(NA, dim = c(n.loc, n.col, len.tau.comp))
    for(i in 1:n.loc){
      for(l in 1:n.col){
        beta.hat[i,l,] <- approx(x = tau.post[k,i,],
                                 y = coefgp$beta.samp[,l,k],
                                 xout = tau.grid)$y
      }
    }
    return(beta.hat)
  }
  
  beta.hat.list <- foreach(k=1:nmc) %dopar% {betahat(k, n.test, p+1, tau.len)}
  stopCluster(cl)
  
  beta.mean <- Reduce('+', beta.hat.list)/nmc
  rm(beta.hat.list)
  gc()
  
  
  cond.qt <- matrix(NA, nrow = n.test, ncol = tau.len)
  x.test.inter <- cbind(1, x.test)
  
  for(i in 1:n.test){
    for(k in 1:tau.len){
      cond.qt[i,k] <- sum(x.test.inter[i,]*beta.mean[i,,k])
    }
  }
  
  jsqr.test.err <- abs(cond.qt - qt.test)
  print(datanum)
}