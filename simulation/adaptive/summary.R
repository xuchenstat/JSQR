library(plyr)
library(quantreg)
library(splines)
library(kernlab)
library(Matrix)


dyn.load('/jsqr/jsqrgp/jsqrgp.so')
source('/jsqr/jsqrgp/utility.R')
source('/jsqr/jsqrgp/jsqrgp.R')


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


alpha.true <- phi.true <- rep(NA, 100)
alpha.qt <- phi.qt <- matrix(NA, nrow = 100, ncol = 4)

cor.true <- matrix(NA, nrow = 100, ncol = 10)
cor.qt <- array(NA, dim = c(100, 10, 4))

for (datanum in 401:500){
  
  set.seed(datanum)    
  
  n <- 500; p <- 1; location <- matrix(runif(2*n), nrow = n)
  distance = as.vector(dist(location, method = 'euclidean'))
  kappa <- 2; cte <- 2^(1-kappa)/gamma(kappa)
  phi.lb <- uniroot(function(u) sapply(u, matern, x = max(distance)/4, kappa = kappa, cte = cte) - 0.05, interval = c(1e-4, 10*max(distance)), tol = 1e-10)$root
  phi.ub <- 3*phi.lb
  
  phi <- runif(1, min = phi.lb, max = phi.ub); 
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
  
  alpha <- runif(1)
  
  z <- t(chol(alpha*C+(1-alpha)*diag(n)))%*%rnorm(n)
  tau <- pnorm(z)
  x <- matrix(runif(n*p,-1,1), nrow = n)
  y <- -3*(tau-0.5)*log(tau*(1-tau))-4*(tau-0.5)^2*log(tau*(1-tau))*x
  

  parmgsamp <- matrix(jsqr.fit$parmgsamp, nrow = 2e3, byrow = TRUE)
  alphasamp <- exp(jsqr.fit$lkappasamp[seq(1002,2000,2)])/exp(parmgsamp[seq(1002,2000,2),15]) 

  phi.grid <- seq(phi.lb, phi.ub, length.out = 10)
  phisamp <- phi.grid[jsqr.fit$phisamp[seq(1002,2000,2)]+1]
  
  alpha.true[datanum - 400] <- alpha
  phi.true[datanum - 400] <- phi
  
  sample.ix <- sample(500,5)
  distance = as.vector(dist(location[sample.ix,], method = 'euclidean'))
  cor.true[datanum - 400,] <- alpha*sapply(distance, matern, phi = phi, kappa = kappa, cte = cte)
  
  
  alpha.qt[datanum - 400,] <- c(quantile(alphasamp, probs = c(0.025,0.5,0.975)), mean(alphasamp))
  phi.qt[datanum - 400,] <- c(quantile(phisamp, probs = c(0.025,0.5,0.975)), mean(phisamp))
  corsamp <- matrix(NA, 500, 10)
  for(i in 1:500){
    corsamp[i,] <- alphasamp[i]*sapply(distance, matern, phi = phisamp[i], kappa = kappa, cte = cte)
  }
  cor.qt[datanum - 400,,] <- t(rbind(apply(corsamp, 2, quantile, probs = c(0.025,0.5,0.975)), apply(corsamp, 2, mean)))
}


length(which(sign((alpha.qt[,1] - alpha.true)*(alpha.qt[,3] - alpha.true))<0))
mean(alpha.qt[,3] - alpha.qt[,1])

length(which(sign((phi.qt[,1] - phi.true)*(phi.qt[,3] - phi.true))<0))
mean(phi.qt[,3] - phi.qt[,1])

length(which(sign((alpha.qt[,1] - alpha.true)*(alpha.qt[,3] - alpha.true))<0))
mean(alpha.qt[,3] - alpha.qt[,1])

mean(abs(alpha.qt[,4] - alpha.true))
sd(abs(alpha.qt[,4] - alpha.true))
length(which(sign((alpha.qt[,1] - alpha.true)*(alpha.qt[,3] - alpha.true))<0))
mean(alpha.qt[,3] - alpha.qt[,1])

mean(abs(phi.qt[,4] - phi.true))
sd(abs(phi.qt[,4] - phi.true))
length(which(sign((phi.qt[,1] - phi.true)*(phi.qt[,3] - phi.true))<0))
mean(phi.qt[,3] - phi.qt[,1])

mae.cor <- ci.len <- matrix(NA, 100, 10)
mae.ci <- rep(NA, 100)

for(i in 1:100){
  mae.cor[i,] <- abs(cor.true[i,] - cor.qt[i,,4])
  mae.ci[i] <- length(which(sign((cor.qt[i,,1] - cor.true[i,])*(cor.qt[i,,3] - cor.true[i,]))<0))
  ci.len[i,] <- cor.qt[i,,3] - cor.qt[i,,1]
}
mean(mae.cor)
sd(mae.cor)
mean(ci.len)
sum(mae.ci)/1000
