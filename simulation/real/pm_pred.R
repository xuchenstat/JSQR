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
library(qrjoint)
library(mvtnorm)
library(doParallel)
library(foreach)

source('/jsqr/jsqrgp/utility.R')

alp.test.loss <- ald.test.loss <- kb.test.loss <- 
  jsqr.test.loss <- jqr.test.loss <- array(NA, dim = c(10, n.test, length(tau.grid)))

alp.train.loss <- ald.train.loss <- kb.train.loss <- 
  jsqr.train.loss <- jqr.train.loss <- array(NA, dim = c(10, n.train, length(tau.grid)))

for(datanum in 1:10){
  
  set.seed(datanum)
  
  train.ix <- sample(n, n.train)
  test.ix <- (1:n)[-train.ix]
  
  y.train <- pm[train.ix,]$pm
  x.train <- pm[train.ix, c('LCYPOP00', 'LDISTA1', 'NLCD_URB', 'LS25_E10', 'ELEV_NED')]
  distance <- as.vector(dist(coords[train.ix,], method = 'euclidean'))
  
  y.test <- pm[test.ix,]$pm
  x.test <- pm[test.ix, c('LCYPOP00', 'LDISTA1', 'NLCD_URB', 'LS25_E10', 'ELEV_NED')]
  
  ######################## jqr
  file.name <- paste('jqrpm', datanum, '.RData', sep = '')
  load(file.name)

  coef.jqr <- coef(jqr.fit, plot = FALSE, nmc = 500)
  coef.jqr.mean <- apply(coef.jqr$beta.samp[as.integer(tau.grid/0.01)+1,,], c(1,2), mean)

  jqr.train.pred <- apply(coef.jqr.mean, 1, function(u) as.matrix(cbind(1,x.train))%*%u)
  jqr.test.pred <- apply(coef.jqr.mean, 1, function(u) as.matrix(cbind(1,x.test))%*%u)

  jqr.train.err <- apply(jqr.train.pred, 2, function(u) y.train - u)
  jqr.test.err <- apply(jqr.test.pred, 2, function(u) y.test - u)
   
  ######################## kb
  file.name <- paste('kbpm', datanum, '.RData', sep = '')
  load(file.name)

  coef.kb <- t(kb.fit$coefficients)

  kb.train.pred <- apply(coef.kb, 1, function(u) as.matrix(cbind(1,x.train))%*%u)
  kb.test.pred <- apply(coef.kb, 1, function(u) as.matrix(cbind(1,x.test))%*%u)

  kb.train.err <- apply(kb.train.pred, 2, function(u) y.train - u)
  kb.test.err <- apply(kb.test.pred, 2, function(u) y.test - u)

  ####################### ald
  file.name <- paste('aldpm', datanum, '.RData', sep = '')
  load(file.name)

  coef.ald <- matrix(NA, nrow = length(tau.grid), ncol = p+1)

  for(i in 1:length(tau.grid)){
    coef.ald[i,] <- apply(ald.fit[[i]]$betadraw[501:1000,], 2, mean)
  }

  ald.train.pred <- apply(coef.ald, 1, function(u) as.matrix(cbind(1,x.train))%*%u)
  ald.test.pred <- apply(coef.ald, 1, function(u) as.matrix(cbind(1,x.test))%*%u)

  ald.train.err <- apply(ald.train.pred, 2, function(u) y.train - u)
  ald.test.err <- apply(ald.test.pred, 2, function(u) y.test - u)

  ######################## alp
  file.name <- paste('alppm', datanum, '.RData', sep = '')
  load(file.name)

  coef.alp <- array(NA, dim = c(500, length(tau.grid), p+1))
  for(taunum in 1:length(tau.grid)){
    coef.alp[,taunum,] <- matrix(alp.fit[[taunum]]$betasamp[((p+1)*500+1):((p+1)*1000)], nrow = 500, ncol = p+1, byrow = TRUE)
  }
  coef.alp.mean <- apply(coef.alp, c(2,3), mean)

  alp.train.pred <- apply(coef.alp.mean, 1, function(u) as.matrix(cbind(1,x.train))%*%u)
  alp.test.pred <- apply(coef.alp.mean, 1, function(u) as.matrix(cbind(1,x.test))%*%u)

  alp.train.err <- apply(alp.train.pred, 2, function(u) y.train - u)
  alp.test.err <- apply(alp.test.pred, 2, function(u) y.test - u)

  ######################## alpadj
  eps <- Z <- W <- array(NA, dim = c(length(tau.grid), 500, n.train))
  quan.adj <- array(NA, dim = c(length(tau.grid), 500, n.test))
  quan.train.adj <- array(NA, dim = c(length(tau.grid), 500, n.train))

  kappa = 2
  cte <- 2^(1-kappa)/gamma(kappa)
  phi.lb <- uniroot(function(u) sapply(u, matern, x = max(distance)/4, kappa = kappa, cte = cte) - 0.05, interval = c(1e-4, 10*max(distance)), tol = 1e-10)$root
  phi.ub <- 3*phi.lb
  phi.grid <- seq(phi.lb, phi.ub, length.out = nphi)
  cov.grid <- lapply(phi.grid, covmat, n = n.train, distance = distance, kernel = 'matern', kappa = kappa, cte = cte)
  Sigma.inv <- lapply(cov.grid, function(u) solve(u+diag(1e-13,n.train)))

  for(taunum in 1:length(tau.grid)){
    tau1 <- (1 - 2*tau.grid[taunum])/tau.grid[taunum]/(1 - tau.grid[taunum])
    tau2 <- 2/tau.grid[taunum]/(1 - tau.grid[taunum])
    tau3 <- sqrt(pi/2/tau.grid[taunum]/(1 - tau.grid[taunum]))
    eps[taunum,,] <- t(apply(coef.alp[,taunum,], 1, function(u) y.train - as.matrix(cbind(1,x.train))%*%u))
    xisamp <- matrix(alp.fit[[taunum]]$xisamp[(n.train*500+1):(n.train*1000)], nrow = 500, ncol = n.train, byrow = TRUE)

    Z[taunum,,] <- (eps[taunum,,] - tau1*xisamp)/sqrt(tau2/alp.fit[[taunum]]$etasamp[(500+1):1000]*xisamp)

    for(iter in 1:500){
      alpha.samp <- alp.fit[[taunum]]$alphasamp[500+iter]
      Omega <- solve(alpha.samp/(1 - alpha.samp)*diag(n.train) + Sigma.inv[[alp.fit[[taunum]]$phisamp[500+iter]+1]])
      W[taunum, iter, ] <- rmvnorm(1, mean = sqrt(alpha.samp)/(1 - alpha.samp)*Omega%*%Z[taunum,iter,], sigma = (Omega + t(Omega))/2)

      quan.train.adj[taunum, iter, ] <- sqrt(alpha.samp)*tau3/alp.fit[[taunum]]$etasamp[500+iter]*W[taunum, iter, ]

      dist.full <- as.vector(dist(coords[c(test.ix,train.ix),], method = 'euclidean'))
      cov.full <- covmat(phi.grid[alp.fit[[taunum]]$phisamp[500+iter]+1], n = n, distance = dist.full, kernel = 'matern', kappa = kappa, cte = cte)

      quan.adj[taunum,iter,] <- sqrt(alpha.samp)*tau3/alp.fit[[taunum]]$etasamp[500+iter]*cov.full[1:n.test, (n.test+1):n]%*%Sigma.inv[[alp.fit[[taunum]]$phisamp[500+iter]+1]]%*%W[taunum, iter, ]

    }
  }

  alpadj.train.pred <- t(apply(quan.train.adj, c(1,3), mean)) + alp.train.pred
  alpadj.train.err <- apply(alpadj.train.pred, 2, function(u) y.train - u)

  alpadj.test.pred <- t(apply(quan.adj, c(1,3), mean)) + alp.test.pred
  alpadj.test.err <- apply(alpadj.test.pred, 2, function(u) y.test - u)
  
  ######################## jsqr
  file.name <- paste('jsqrpm', datanum, '.RData', sep = '')
  load(file.name)

  nsamp <- jsqr.fit$dim[11]; burn.perc <- 0.5; nmc <- 5e2; m <- jsqr.fit$dim[5]
  usamp <- matrix(jsqr.fit$usamp,ncol=length(train.ix),byrow=T)
  parmgsamp <- matrix(jsqr.fit$parmgsamp, nrow = nsamp, byrow = TRUE)
  coefgp <- coef.qrjointgp(jsqr.fit, plot=FALSE, nmc = nmc, burn.perc = burn.perc, reduce = FALSE)
  ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 1)[-1]))

  nphi <- 10; kappa <- 2
  cte <- 2^(1-kappa)/gamma(kappa)
  phi.lb <- uniroot(function(u) sapply(u, matern, x = max(distance)/4, kappa = kappa, cte = cte) - 0.05, interval = c(1e-4, 10*max(distance)), tol = 1e-10)$root
  phi.ub <- 3*phi.lb
  phi.grid <- seq(phi.lb, phi.ub, length.out = nphi)

  tau.g <- jsqr.fit$tau.g
  z.g <- qnorm(tau.g)
  len.reg.ix <- length(z.g)
  tau.post <- array(NA, dim = c(nmc, length(test.ix), len.reg.ix))
  tau.train.post <- array(NA, dim = c(nmc, length(train.ix), len.reg.ix))

  cov.grid <- array(NA, dim = c(nphi, length(train.ix), length(train.ix)))
  for(i in 1:nphi){
    cov.grid[i,,] <- covmat(phi = phi.grid[i], distance = distance, n = length(train.ix), kernel = 2, kappa = kappa, cte = cte)
  }
  distance.test <- t(apply(coords[test.ix,], 1, function(u) apply(coords[train.ix,], 1, function(v) sqrt(sum((u-v)^2)))))

  for(i in 1:nmc){
    alpha <- exp(jsqr.fit$lkappasamp[ss[i]])/exp(parmgsamp[ss[i],(m+1)*(p+1)+1])
    sigma22 <- alpha*cov.grid[jsqr.fit$phisamp[ss[i]]+1,,]+(1-alpha)*diag(n.train)
    sigma22.inv <- chol2inv(chol(sigma22))
    zsamp <- qnorm(usamp[ss[i],])
    for(j in 1:n.test){
      sigma11 <- 1
      sigma12 <- alpha*corfun(phi = phi.grid[jsqr.fit$phisamp[ss[i]]+1], distance = distance.test[j,], kernel = 2, kappa = kappa, cte = cte)
      sigma12.22.inv <- sigma12%*%sigma22.inv
      cond.var <- as.numeric(sigma11 - tcrossprod(sigma12.22.inv, sigma12))
      cond.mean <- c(sigma12.22.inv %*% zsamp)
      tau.post[i,j,] <- pnorm((z.g - cond.mean)/sqrt(cond.var))
    }

    for(j in 1:n.train){
      sigma11 <- 1
      sigma12 <- sigma22[j, 1:n.train, drop = FALSE]
      sigma12[j] <- alpha
      sigma12.22.inv <- sigma12%*%sigma22.inv
      cond.var <- as.numeric(sigma11 - tcrossprod(sigma12.22.inv, sigma12))
      cond.mean <- c(sigma12.22.inv %*% zsamp)
      tau.train.post[i,j,] <- pnorm((z.g - cond.mean)/sqrt(cond.var))
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

  jsqr.test.err <- apply(cond.qt, 2, function(u) y.test - u)

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

  jsqr.train.err <- apply(cond.train.qt, 2, function(u) y.train - u)
  
  for(i in 1:length(tau.grid)){
    jqr.train.loss[datanum,,i] <- check.loss(jqr.train.err[,i], tau.grid[i])
    jqr.test.loss[datanum,,i] <- check.loss(jqr.test.err[,i], tau.grid[i])
    
    jsqr.train.loss[datanum,,i] <- check.loss(jsqr.train.err[,i], tau.grid[i])
    jsqr.test.loss[datanum,,i] <- check.loss(jsqr.test.err[,i], tau.grid[i])

    kb.train.loss[datanum,,i] <- check.loss(kb.train.err[,i], tau.grid[i])
    kb.test.loss[datanum,,i] <- check.loss(kb.test.err[,i], tau.grid[i])
    

    alp.train.loss[datanum,,i] <- check.loss(alpadj.train.err[,i], tau.grid[i])
    alp.test.loss[datanum,,i] <- check.loss(alpadj.test.err[,i], tau.grid[i])
    
    jsqr.train.loss[datanum,,i] <- check.loss(jsqr.train.err[,i], tau.grid[i])
    jsqr.test.loss[datanum,,i] <- check.loss(jsqr.test.err[,i], tau.grid[i])
  }
  print(datanum)
}

# save(jqr.train.loss, file = 'jqrtrainloss.RData')
# save(jqr.test.loss, file = 'jqrtestloss.RData')
# save(kb.train.loss, file = 'kbtrainloss.RData')
# save(kb.test.loss, file = 'kbtestloss.RData')
# save(ald.train.loss, file = 'aldtrainloss.RData')
# save(ald.test.loss, file = 'aldtestloss.RData')
# save(alp.train.loss, file = 'alptrainloss.RData')
# save(alp.test.loss, file = 'alptestloss.RData')
# save(jsqr.train.loss, file = 'jsqrtrainloss.RData')
# save(jsqr.test.loss, file = 'jsqrtestloss.RData')
