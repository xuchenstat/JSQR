library(ald)
library(plyr)
library(quantreg)
library(splines)
library(kernlab)
library(Matrix)
library(parallel)

dyn.load('/lum2012/alpnew.so')
source('/lum2012/alQR.R')

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

for (datanum in 1:100){

set.seed(datanum)    

p.ald <- 0.4
tau.ald <- 1

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

z.ald <- t(chol(alpha*C+(1-alpha)*diag(n)))%*%rnorm(n)
xi.ald <- rexp(n, tau.ald)
tau <- pALD(sqrt(2*xi.ald/tau.ald/p.ald/(1-p.ald))*z.ald+(1-2*p.ald)/p.ald/(1-p.ald)*xi.ald, mu = 0, sigma = 1/tau.ald)


x <- matrix(runif(n*p,-1,1), nrow = n)
y <- -3*(tau-0.5)*log(tau*(1-tau))-4*(tau-0.5)^2*log(tau*(1-tau))*x

tau.grid <- c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99)


alp.fit <- mclapply(tau.grid, alqr, x = cbind(1,x), y = y, nphi = 10, alpha = 0.5, update.alpha = TRUE, sdlalpha = 2, decaypar = c(0, max(distance)),
                  distance = distance, kappa = 2, kernel = 'matern', nsamp = 1e3, thin = 10, verbose = TRUE, mc.preschedule = TRUE, mc.cores = 5)

for(taunum in 1:length(tau.grid)){
  alp.fit[[taunum]]$S <- NULL
  alp.fit[[taunum]]$x <- NULL
  alp.fit[[taunum]]$y <- NULL
  
  alp.fit[[taunum]]$betasamp <- alp.fit[[taunum]]$betasamp[((p+1)*500+1):((p+1)*1000)]
  alp.fit[[taunum]]$xisamp <- alp.fit[[taunum]]$xisamp[(n*500+1):(n*1000)]
  alp.fit[[taunum]]$etasamp <- alp.fit[[taunum]]$etasamp[(500+1):1000]
  alp.fit[[taunum]]$phisamp <- alp.fit[[taunum]]$phisamp[(500+1):1000]
  alp.fit[[taunum]]$alphasamp <- alp.fit[[taunum]]$alphasamp[(500+1):1000]
}

save(alp.fit, file = paste('alpfit', datanum, '.RData', sep = ''))
}
