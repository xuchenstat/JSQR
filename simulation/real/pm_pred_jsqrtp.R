pm <- read.csv('pm.csv') 

n <- nrow(pm)
n.train <- 270
n.test <- n - n.train
tau.grid <- c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99)
coords <- pm[,c('x','y')]
p <- 5

check.loss <- function(r, tau){
  return(r*(tau - (r<0)))
}

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

nphi <- 10

library(kernlab)
library(Matrix)
library(mvtnorm)
library(doParallel)
library(foreach)
source('/jsqr/jsqrtp/utility.R')

jsqrtp.test.loss <- array(NA, dim = c(10, n.test, length(tau.grid)))
jsqrtp.train.loss <- array(NA, dim = c(10, n.train, length(tau.grid)))

for(datanum in 1:10){
  
  set.seed(datanum)
  
  train.ix <- sample(n, n.train)
  test.ix <- (1:n)[-train.ix]
  
  y.train <- pm[train.ix,]$pm
  x.train <- pm[train.ix, c('LCYPOP00', 'LDISTA1', 'NLCD_URB', 'LS25_E10', 'ELEV_NED')]
  distance <- as.vector(dist(coords[train.ix,], method = 'euclidean'))
  
  y.test <- pm[test.ix,]$pm
  x.test <- pm[test.ix, c('LCYPOP00', 'LDISTA1', 'NLCD_URB', 'LS25_E10', 'ELEV_NED')]
  
  ######################## jsqrtp
  file.name <- paste('jsqrtppm', datanum, '.RData', sep = '')

  load(file.name)
  
  nsamp <- jsqr.fit$dim[11]; burn.perc <- 0.5; nmc <- 5e2; m <- jsqr.fit$dim[5]
  usamp <- matrix(jsqr.fit$usamp,ncol=length(train.ix),byrow=T)
  parmgsamp <- matrix(jsqr.fit$parmgsamp, nrow = nsamp, byrow = TRUE)
  coefgp <- coef.qrjointtp(jsqr.fit, plot=FALSE, nmc = nmc, burn.perc = burn.perc, reduce = FALSE)          
  ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 1)[-1]))
  
  nphi <- 10; kappa <- 2
  cte <- 2^(1-kappa)/gamma(kappa)
  phi.lb <- uniroot(function(u) sapply(u, matern, x = max(distance)/4, kappa = kappa, cte = cte) - 0.05, interval = c(1e-4, 10*max(distance)), tol = 1e-10)$root
  phi.ub <- 3*phi.lb
  phi.grid <- seq(phi.lb, phi.ub, length.out = nphi) 
  
  tau.g <- jsqr.fit$tau.g
  len.reg.ix <- length(tau.g)
  tau.post <- array(NA, dim = c(nmc, length(test.ix), len.reg.ix))
  tau.train.post <- array(NA, dim = c(nmc, length(train.ix), len.reg.ix))
  
  cov.grid <- array(NA, dim = c(nphi, length(train.ix), length(train.ix)))
  for(i in 1:nphi){
    cov.grid[i,,] <- covmat(phi = phi.grid[i], distance = distance, n = length(train.ix), kernel = 2, kappa = kappa, cte = cte)
  }
  distance.test <- t(apply(coords[test.ix,], 1, function(u) apply(coords[train.ix,], 1, function(v) sqrt(sum((u-v)^2)))))
  
  for(i in 1:nmc){
    psi <- plogis(jsqr.fit$lpsisamp[ss[i]])*18+2
    z.g <- qt(tau.g, df = psi)
    zsamp <- qt(usamp[ss[i],], df = psi)
    
    alpha <- exp(jsqr.fit$lkappasamp[ss[i]])/exp(parmgsamp[ss[i],(m+1)*(p+1)+1])
    sigma22 <- alpha*cov.grid[jsqr.fit$phisamp[ss[i]]+1,,]+(1-alpha)*diag(n.train)
    sigma22.inv <- chol2inv(chol(sigma22))
    dmaha <- t(zsamp)%*%sigma22.inv%*%zsamp
    
    for(j in 1:n.test){
      sigma11 <- 1
      sigma12 <- alpha*corfun(phi = phi.grid[jsqr.fit$phisamp[ss[i]]+1], distance = distance.test[j,], kernel = 2, kappa = kappa, cte = cte)
      sigma12.22.inv <- sigma12%*%sigma22.inv
      cond.var <- as.numeric((psi + dmaha)*(sigma11 - tcrossprod(sigma12.22.inv, sigma12))/(psi + n.train))
      cond.mean <- c(sigma12.22.inv %*% zsamp)
      tau.post[i,j,] <- pt((z.g - cond.mean)/sqrt(cond.var), df = psi + n.train)
    }
    
    for(j in 1:n.train){
      sigma11 <- 1
      sigma12 <- sigma22[j, 1:n.train, drop = FALSE]
      sigma12[j] <- alpha
      sigma12.22.inv <- sigma12%*%sigma22.inv
      cond.var <- as.numeric((psi + dmaha)*(sigma11 - tcrossprod(sigma12.22.inv, sigma12))/(psi + n.train))
      cond.mean <- c(sigma12.22.inv %*% zsamp)
      tau.train.post[i,j,] <- pt((z.g - cond.mean)/sqrt(cond.var), df = psi + n.train)
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
  
  beta.hat.list <- foreach(k=1:nmc) %dopar% {betahat(k, length(test.ix), p+1, length(tau.grid))}
  stopCluster(cl)
  
  beta.mean <- Reduce('+', beta.hat.list)/nmc
  rm(beta.hat.list)
  gc()
  
  
  cond.qt <- matrix(NA, nrow = length(test.ix), ncol = length(tau.grid))
  x.test.inter <- cbind(1, x.test)
  
  for(i in 1:length(test.ix)){
    for(k in 1:length(tau.grid)){
      cond.qt[i,k] <- sum(x.test.inter[i,]*beta.mean[i,,k])
    }
  }
  
  jsqrtp.test.err <- apply(cond.qt, 2, function(u) y.test - u)
  
  cl <- makeCluster(4)
  registerDoParallel(cl)

  betahat <- function(k, n.loc, n.col, len.tau.comp){
    beta.hat <- array(NA, dim = c(n.loc, n.col, len.tau.comp))
    for(i in 1:n.loc){
      for(l in 1:n.col){
        beta.hat[i,l,] <- approx(x = tau.train.post[k,i,],
                                 y = coefgp$beta.samp[,l,k],
                                 xout = tau.grid)$y
      }
    }
    return(beta.hat)
  }

  beta.hat.list <- foreach(k=1:nmc) %dopar% {betahat(k, n.train, p+1, length(tau.grid))}
  stopCluster(cl)

  beta.mean <- Reduce('+', beta.hat.list)/nmc
  rm(beta.hat.list)
  gc()


  cond.train.qt <- matrix(NA, nrow = n.train, ncol = length(tau.grid))
  x.train.inter <- cbind(1, x.train)

  for(i in 1:n.train){
    for(k in 1:length(tau.grid)){
      cond.train.qt[i,k] <- sum(x.train.inter[i,]*beta.mean[i,,k])
    }
  }

  jsqrtp.train.err <- apply(cond.train.qt, 2, function(u) y.train - u)
  
  for(i in 1:length(tau.grid)){
    jsqrtp.train.loss[datanum,,i] <- check.loss(jsqrtp.train.err[,i], tau.grid[i])
    jsqrtp.test.loss[datanum,,i] <- check.loss(jsqrtp.test.err[,i], tau.grid[i])
  }
  
  print(datanum)
}