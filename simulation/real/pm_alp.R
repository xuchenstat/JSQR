pm <- read.csv('pm.csv') 
n <- nrow(pm)
n.train <- 270
n.test <- n - n.train

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

tau.grid <- c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99)

for(datanum in 1:10){
  
  set.seed(datanum)
  
  train.ix <- sample(n, n.train)
  test.ix <- (1:n)[-train.ix]
  
  y.train <- pm[train.ix,]$pm
  x.train <- pm[train.ix, c('LCYPOP00', 'LDISTA1', 'NLCD_URB', 'LS25_E10', 'ELEV_NED')]
  coords <- pm[,c('x','y')]
  
  distance <- as.vector(dist(coords[train.ix,], method = 'euclidean'))
  
  alp.fit <- mclapply(tau.grid, alqr, x = as.matrix(cbind(1, x.train)), y = y.train, nphi = 10, alpha = 0.5, update.alpha = TRUE, sdlalpha = 2, decaypar = c(0, max(distance)),
                    distance = distance, kappa = 2, kernel = 'matern', nsamp = 1e3, thin = 10, verbose = TRUE, mc.preschedule = TRUE, mc.cores = 13)
  
  save(alp.fit, file = paste('alppm', datanum, '.RData', sep = ''))
}