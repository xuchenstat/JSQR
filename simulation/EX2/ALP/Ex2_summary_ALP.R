tau.grid <- c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99)

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

p <- length(b.half)

btrue <- matrix(NA, nrow = length(tau.grid), ncol = p+1)

btrue[,1] <- sapply(tau.grid, function(u) integrate(b0.prime, lower = 0.5, upper = u, rel.tol = 1e-10)$value + b0.half)
btrue[,2] <- sapply(tau.grid, function(u) integrate(b.prime, lower = 0.5, upper = u, j = 1, rel.tol = 1e-10)$value + b.half[1])
btrue[,3] <- sapply(tau.grid, function(u) integrate(b.prime, lower = 0.5, upper = u, j = 2, rel.tol = 1e-10)$value + b.half[2])
btrue[,4] <- sapply(tau.grid, function(u) integrate(b.prime, lower = 0.5, upper = u, j = 3, rel.tol = 1e-10)$value + b.half[3])
btrue[,5] <- sapply(tau.grid, function(u) integrate(b.prime, lower = 0.5, upper = u, j = 4, rel.tol = 1e-10)$value + b.half[4])
btrue[,6] <- sapply(tau.grid, function(u) integrate(b.prime, lower = 0.5, upper = u, j = 5, rel.tol = 1e-10)$value + b.half[5])
btrue[,7] <- sapply(tau.grid, function(u) integrate(b.prime, lower = 0.5, upper = u, j = 6, rel.tol = 1e-10)$value + b.half[6])
btrue[,8] <- sapply(tau.grid, function(u) integrate(b.prime, lower = 0.5, upper = u, j = 7, rel.tol = 1e-10)$value + b.half[7])



source('/jsqr/jsqrgp/utility.R')

datanums <- 1:100
mean.abs.error <- array(NA, dim = c(5, 100, length(tau.grid), p+1))
cp <- array(NA, dim = c(5, 100, length(tau.grid), p+1))

for(datanum in datanums){
  filename <- paste('jsqrfit', datanum + 200, '.RData', sep = '')
  load(filename)
  
  coef.jsqr <- coef.qrjointgp(jsqr.fit, nmc = 500, plot = FALSE)
  
  post.mean <- apply(coef.jsqr$beta.samp[as.integer(tau.grid/0.01)+1,,], c(1,2), mean)
  post.lq <- apply(coef.jsqr$beta.samp[as.integer(tau.grid/0.01)+1,,], c(1,2), quantile, 0.025)
  post.uq <- apply(coef.jsqr$beta.samp[as.integer(tau.grid/0.01)+1,,], c(1,2), quantile, 0.975)
  
  mean.abs.error[1, datanum,,] <- abs(btrue - post.mean)
  cp[1, datanum,,] <- (btrue > post.lq) & (btrue < post.uq)
}


library(qrjoint)
for(datanum in datanums){
  filename <- paste('jqrfit', datanum + 200, '.RData', sep = '')
  load(filename)
  
  coef.jqr <- coef(jqr.fit, nmc = 500, plot = FALSE)  
  
  post.mean <- apply(coef.jqr$beta.samp[as.integer(tau.grid/0.01)+1,,], c(1,2), mean)
  post.lq <- apply(coef.jqr$beta.samp[as.integer(tau.grid/0.01)+1,,], c(1,2), quantile, 0.025)
  post.uq <- apply(coef.jqr$beta.samp[as.integer(tau.grid/0.01)+1,,], c(1,2), quantile, 0.975)
  
  mean.abs.error[2, datanum,,] <- abs(btrue - post.mean)
  cp[2, datanum,,] <- (btrue > post.lq) & (btrue < post.uq)
}


for(datanum in datanums){
  coef.lq <- coef.uq <- matrix(NA, nrow = length(tau.grid), ncol = p+1)
  filename <- paste('kbfit', datanum + 200, '.RData', sep = '')
  load(filename)
  
  coef.kb <- t(kb.fit$coefficients)
  kb.summary <- summary(kb.fit)
  
  for(i in 1:length(tau.grid)){
    coef.lq[i,] <- kb.summary[[i]]$coefficients[,2]
    coef.uq[i,] <- kb.summary[[i]]$coefficients[,3]
  }
  
  mean.abs.error[3, datanum,,] <- abs(btrue - coef.kb)
  cp[3, datanum,,] <- (btrue > coef.lq) & (btrue < coef.uq)
}


for(datanum in datanums){
  coef.lq <- coef.uq <- matrix(NA, nrow = length(tau.grid), ncol = p+1)
  filename <- paste('alpfit', datanum + 200, '.RData', sep = '')
  load(filename)
  
  coef.alp <- array(NA, dim = c(500, length(tau.grid), p+1))
  for(taunum in 1:length(tau.grid)){
    coef.alp[,taunum,] <- matrix(alp.fit[[taunum]]$betasamp, nrow = 500, ncol = p + 1, byrow = TRUE)
  }
  coef.alp.mean <- apply(coef.alp, c(2,3), mean)
  
  coef.lq <- apply(coef.alp, c(2,3), quantile, 0.025)
  coef.uq <- apply(coef.alp, c(2,3), quantile, 0.975)
  
  mean.abs.error[4, datanum,,] <- abs(btrue - coef.alp.mean)
  cp[4, datanum,,] <- (btrue > coef.lq) & (btrue < coef.uq)
}


for(datanum in datanums){
  coef.ald <- coef.lq <- coef.uq <- matrix(NA, nrow = length(tau.grid), ncol = p+1)
  filename <- paste('aldfit', datanum + 200, '.RData', sep = '')
  load(filename)
  
  for(i in 1:length(tau.grid)){
    coef.ald[i,] <- apply(ald.fit[[i]]$betadraw[501:1000,], 2, mean)
    coef.lq[i,] <- apply(ald.fit[[i]]$betadraw[501:1000,], 2, quantile, probs = 0.025)
    coef.uq[i,] <- apply(ald.fit[[i]]$betadraw[501:1000,], 2, quantile, probs = 0.975)
  }
  
  mean.abs.error[5, datanum,,] <- abs(btrue - coef.ald)
  cp[5, datanum,,] <- (btrue > coef.lq) & (btrue < coef.uq)
}


save(mean.abs.error, file = 'ex2alpmae.RData')
save(cp, file = 'ex2alpcp.RData')


mae.jsqr <- apply(mean.abs.error[1,,,], c(2,3), mean)
mae.jqr <- apply(mean.abs.error[2,,,], c(2,3), mean)
mae.kb <- apply(mean.abs.error[3,,,], c(2,3), mean)
mae.alp <- apply(mean.abs.error[4,,,], c(2,3), mean)


cp.jsqr <- apply(cp[1,,,], c(2,3), mean)
cp.jqr <- apply(cp[2,,,], c(2,3), mean)
cp.kb <- apply(cp[3,,,], c(2,3), mean)
cp.alp <- apply(cp[4,,,], c(2,3), mean)