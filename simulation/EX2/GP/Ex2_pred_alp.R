tau.grid <- c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99)
tau.len <- length(tau.grid)
z.grid <- qnorm(tau.grid)

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

nphi <- 10

library(kernlab)
library(Matrix)
library(mvtnorm)

for(datanum in 301:400){
  
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
  
  location.test <- matrix(runif(2*n.test), nrow = n.test)
  
  z2 <- matrix(rnorm(n.test*p),n.test,p)
  zn2 <- sweep(z2,1,sqrt(apply(z2^2,1,sum)),'/')
  u2 <- rbeta(n.test,p,1)
  x.test <- sweep(zn2,1,u2,'*')
  distance.test <- t(apply(location.test, 1, function(u) apply(location, 1, function(v) sqrt(sum((u-v)^2)))))
  
  sigma22 <- alpha*C+(1-alpha)*diag(n)
  sigma22.inv <- chol2inv(chol(sigma22))
  
  tau.test <- qt.test <- matrix(NA, nrow = n.test, ncol = tau.len)
  
  for(i in 1:n.test){
    sigma11 <- 1
    sigma12 <- alpha*corfun(phi = phi, distance = distance.test[i,], kernel = 2, kappa = kappa, cte = cte)
    sigma12.22.inv <- sigma12%*%sigma22.inv
    cond.var <- as.numeric(sigma11 - tcrossprod(sigma12.22.inv, sigma12))
    cond.mean <- c(sigma12.22.inv %*% z)
    tau.test[i,] <- pnorm(cond.mean + sqrt(cond.var) * z.grid)
    for(k in 1:tau.len){
      b0 <- integrate(b0.prime, 0.5, tau.test[i,k], rel.tol = 1e-10)$value + b0.half
      b <- sapply(1:p, function(u) integrate(b.prime, lower = 0.5, upper = tau.test[i,k], j = u, rel.tol = 1e-10)$value) + b.half
      qt.test[i,k] <- b0 + b%*%x.test[i,]
    }
  }
  
  
  file.name <- paste('alpfit', datanum, '.RData', sep = '')
  load(file.name)
  
  coef.alp <- array(NA, dim = c(500, tau.len, p+1))
  for(taunum in 1:tau.len){
    coef.alp[,taunum,] <- matrix(alp.fit[[taunum]]$betasamp, nrow = 500, ncol = p+1, byrow = TRUE)
  }
  coef.alp.mean <- apply(coef.alp, c(2,3), mean)
  
  alp.test.pred <- apply(coef.alp.mean, 1, function(u) as.matrix(cbind(1,x.test))%*%u)
  
  eps <- Z <- W <- array(NA, dim = c(tau.len, 500, n))
  quan.adj <- array(NA, dim = c(tau.len, 500, n.test))
  
  phi.lb <- uniroot(function(u) sapply(u, matern, x = max(distance)/4, kappa = kappa, cte = cte) - 0.05, interval = c(1e-4, 10*max(distance)), tol = 1e-10)$root
  phi.ub <- 3*phi.lb
  phi.grid <- seq(phi.lb, phi.ub, length.out = nphi)
  cov.grid <- lapply(phi.grid, covmat, n = n, distance = distance, kernel = 'matern', kappa = kappa, cte = cte)
  Sigma.inv <- lapply(cov.grid, function(u) solve(u))

  for(taunum in 1:tau.len){
    tau1 <- (1 - 2*tau.grid[taunum])/tau.grid[taunum]/(1 - tau.grid[taunum])
    tau2 <- 2/tau.grid[taunum]/(1 - tau.grid[taunum])
    tau3 <- sqrt(pi/2/tau.grid[taunum]/(1 - tau.grid[taunum]))
    eps[taunum,,] <- t(apply(coef.alp[,taunum,], 1, function(u) y - as.matrix(cbind(1,x))%*%u))
    xisamp <- matrix(alp.fit[[taunum]]$xisamp, nrow = 500, ncol = n, byrow = TRUE)
    
    Z[taunum,,] <- (eps[taunum,,] - tau1*xisamp)/sqrt(tau2/alp.fit[[taunum]]$etasamp*xisamp)
    
    for(iter in 1:500){
      alpha.samp <- alp.fit[[taunum]]$alphasamp[iter]
      Omega <- solve(alpha.samp/(1 - alpha.samp)*diag(n) + Sigma.inv[[alp.fit[[taunum]]$phisamp[iter]+1]])
      W[taunum, iter, ] <- rmvnorm(1, mean = sqrt(alpha.samp)/(1 - alpha.samp)*Omega%*%Z[taunum,iter,], sigma = (Omega + t(Omega))/2)
      
      dist.full <- as.vector(dist(rbind(location.test,location), method = 'euclidean'))
      cov.full <- covmat(phi.grid[alp.fit[[taunum]]$phisamp[iter]+1], n = n + n.test, distance = dist.full, kernel = 'matern', kappa = kappa, cte = cte)
      
      quan.adj[taunum,iter,] <- sqrt(alpha.samp)*tau3/alp.fit[[taunum]]$etasamp[iter]*cov.full[1:n.test, (n.test+1):(n + n.test)]%*%Sigma.inv[[alp.fit[[taunum]]$phisamp[iter]+1]]%*%W[taunum, iter, ]
      if(iter %% 10 == 0) {print(c(taunum, iter))}
    }
  }
  
  alpadj.test.pred <- t(apply(quan.adj, c(1,3), mean)) + alp.test.pred
  alpadj.test.err <- abs(alpadj.test.pred - qt.test)
  
  print(datanum)
}