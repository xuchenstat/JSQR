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

nphi <- 10

library(kernlab)
library(Matrix)
library(mvtnorm)

for(datanum in 101:200){
  
  set.seed(datanum)    
  
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
  
  z <- t(chol(alpha*C+(1-alpha)*diag(n)))%*%rnorm(n)
  tau <- pnorm(z)
  x <- matrix(runif(n*p,-1,1), nrow = n)
  y <- -3*(tau-0.5)*log(tau*(1-tau))-4*(tau-0.5)^2*log(tau*(1-tau))*x
  
  location.test <- matrix(runif(2*n.test), nrow = n.test)
  x.test <- matrix(runif(n.test*p,-1,1), nrow = n.test)
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
    for(j in 1:tau.len){
      qt.test[i,j] <- -3*(tau.test[i,j]-0.5)*log(tau.test[i,j]*(1-tau.test[i,j]))-4*(tau.test[i,j]-0.5)^2*log(tau.test[i,j]*(1-tau.test[i,j]))*x.test[i,]
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