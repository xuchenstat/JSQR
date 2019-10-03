pm <- read.csv('pm.csv') 
n <- nrow(pm)
n.train <- 270
n.test <- n - n.train

library(kernlab)
library(Matrix)
library(quantreg)

dyn.load('/jsqr/jsqrtp/jsqrtp.so')
source('/jsqr/jsqrtp/utility.R')
source('/jsqr/jsqrtp/jsqrtp.R')


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


for(datanum in 1:10){
  
  set.seed(datanum)
  
  train.ix <- sample(n, n.train)
  test.ix <- (1:n)[-train.ix]
  
  y.train <- pm[train.ix,]$pm
  x.train <- pm[train.ix, c('LCYPOP00', 'LDISTA1', 'NLCD_URB', 'LS25_E10', 'ELEV_NED')]
  coords <- pm[,c('x','y')]
  
  distance <- as.vector(dist(coords[train.ix,], method = 'euclidean'))
  
  jsqr.fit <- qrjointtp(x = x.train, y = y.train, par = 'prior', distance = distance,
                        nphi = 10, kernel = 'matern',
                        acpt.target = 0.15, acpt.hptarget = 0.15, kappa = 2, nsamp = 2e3,
                        thin = 10, fbase = 't')
  
  save(jsqr.fit, file = paste('jsqrtppm', datanum, '.RData', sep = ''))
}