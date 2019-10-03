tau.grid <- c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99)

p <- 1

source('/jsqr/jsqrgp/utility.R')

tau.intercept <- function(tau) {return(-3*(tau-0.5)*log(tau*(1-tau)))}
tau.slope <- function(tau) {return(-4*(tau-0.5)^2*log(tau*(1-tau)))}


btrue <- matrix(NA, nrow = length(tau.grid), ncol = p+1)

btrue[,1] <- tau.intercept(tau.grid)
btrue[,2] <- tau.slope(tau.grid)


datanums <- 1:100
mean.abs.error <- array(NA, dim = c(5, 100, length(tau.grid), p+1))
cp <- array(NA, dim = c(5, 100, length(tau.grid), p+1))


for(datanum in datanums){
  filename <- paste('jsqrfit', datanum, '.RData', sep = '')
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
  filename <- paste('jqrfit', datanum, '.RData', sep = '')
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
  filename <- paste('kbfit', datanum, '.RData', sep = '')
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
  filename <- paste('alpfit', datanum, '.RData', sep = '')
  load(filename)
  
  coef.alp <- array(NA, dim = c(500, length(tau.grid), p+1))
  for(taunum in 1:length(tau.grid)){
    coef.alp[,taunum,] <- matrix(alp.fit[[taunum]]$betasamp, nrow = 500, ncol = 2, byrow = TRUE)
  }
  coef.alp.mean <- apply(coef.alp, c(2,3), mean)
  
  coef.lq <- apply(coef.alp, c(2,3), quantile, 0.025)
  coef.uq <- apply(coef.alp, c(2,3), quantile, 0.975)
  
  mean.abs.error[4, datanum,,] <- abs(btrue - coef.alp.mean)
  cp[4, datanum,,] <- (btrue > coef.lq) & (btrue < coef.uq)
}


for(datanum in datanums){
  coef.ald <- coef.lq <- coef.uq <- matrix(NA, nrow = length(tau.grid), ncol = p+1)
  filename <- paste('aldfit', datanum, '.RData', sep = '')
  load(filename)
  
  for(i in 1:length(tau.grid)){
    coef.ald[i,] <- apply(ald.fit[[i]]$betadraw[501:1000,], 2, mean)
    coef.lq[i,] <- apply(ald.fit[[i]]$betadraw[501:1000,], 2, quantile, probs = 0.025)
    coef.uq[i,] <- apply(ald.fit[[i]]$betadraw[501:1000,], 2, quantile, probs = 0.975)
  }
  
  mean.abs.error[5, datanum,,] <- abs(btrue - coef.ald)
  cp[5, datanum,,] <- (btrue > coef.lq) & (btrue < coef.uq)
}


save(mean.abs.error, file = 'ex1alpmae.RData')
save(cp, file = 'ex1alpcp.RData')


mae.jsqr <- apply(mean.abs.error[1,,,], c(2,3), mean)
mae.jqr <- apply(mean.abs.error[2,,,], c(2,3), mean)
mae.kb <- apply(mean.abs.error[3,,,], c(2,3), mean)
mae.alp <- apply(mean.abs.error[4,,,], c(2,3), mean)


cp.jsqr <- apply(cp[1,,,], c(2,3), mean)
cp.jqr <- apply(cp[2,,,], c(2,3), mean)
cp.kb <- apply(cp[3,,,], c(2,3), mean)
cp.alp <- apply(cp[4,,,], c(2,3), mean)