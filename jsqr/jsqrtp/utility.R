library(mvtnorm)
library(plyr)
library(quantreg)
library(splines)
library(kernlab)
library(Matrix)
library(doParallel)
library(foreach)


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


coef.qrjointtp <- function(object, burn.perc = 0.5, nmc = 200, plot = FALSE, show.intercept = TRUE, reduce = TRUE, ...){
  nsamp <- object$dim[11]                      # number of iterations retains in sample
  pars <- matrix(object$parmgsamp, ncol = nsamp) # reorder parameters into npar x nsamp matrix
  ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 1)[-1])) # indices over which to summarizing (exclude burn; keep nmc samples)
  
  n <- object$dim[1]; p <- object$dim[2]; L <- object$dim[3]; mid <- object$dim[4] + 1; nknots <- object$dim[5]; ngrid <- object$dim[6]
  a.sig <- object$hyper[1:2]; a.kap <- matrix(object$hyper[-c(1:2)], nrow = 3)
  tau.g <- object$tau.g; reg.ix <- object$reg.ix
  x.ce <- outer(rep(1, L), attr(object$x, "scaled:center")); x.sc <- outer(rep(1,L), attr(object$x, "scaled:scale"))
  
  base.bundle <- list()
  if(object$fbase.choice == 1){
    base.bundle$q0 <- function(u, nu = Inf) return(1 / (dt(qt(unitFn(u), df = nu), df = nu) * qt(.9, df = nu)))
    base.bundle$Q0 <- function(u, nu = Inf) return(qt(unitFn(u), df = nu) / qt(.9, df = nu))
    base.bundle$F0 <- function(x, nu = Inf) return(pt(x*qt(.9, df = nu), df = nu))
  } else if(object$fbase.choice == 2){
    base.bundle$q0 <- function(u, nu = Inf) return(1 / dlogis(qlogis(unitFn(u))))
    base.bundle$Q0 <- function(u, nu = Inf) return(qlogis(unitFn(u)))
    base.bundle$F0 <- function(x, nu = Inf) return(plogis(x))
  } else {
    base.bundle$q0 <- function(u, nu = Inf) return(1 / (dunif(qunif(u, -1,1), -1,1)))
    base.bundle$Q0 <- function(u, nu = Inf) return(qunif(u, -1,1))
    base.bundle$F0 <- function(x, nu = Inf) return(punif(x, -1,1))
  }
  
  # use estFn to turn posterior on parameters into posterior betas;
  # return 3-D array { L x (p+1) x nsim }
  beta.samp <- sapply(ss, function(p1) estFn(pars[,p1], object$x, object$y, object$gridmats, L, mid, nknots, ngrid, a.kap, a.sig, tau.g, reg.ix, reduce, x.ce, x.sc, base.bundle), simplify="array")
  
  if(reduce) tau.g <- tau.g[reg.ix]
  L <- length(tau.g)
  if(plot){
    nc <- ceiling(sqrt(p+show.intercept)); nr <- ceiling((p+show.intercept)/nc)
    cur.par <- par(no.readonly = TRUE)
    par(mfrow = c(nr, nc))
  }
  
  # store summaries in 3-D array {L x (p+1) x 3}
  beta.hat <- array(0,c(L,p+1,3))
  plot.titles <- c("Intercept", paste('X', 1:p, sep=''))
  beta.hat[,1,] <- getBands(beta.samp[,1,], plot = (plot & show.intercept), add = FALSE, x = tau.g, xlab = expression(tau), ylab = "Coefficient", bty = "n")
  if(plot & show.intercept) {title(main = plot.titles[1])}
  for(j in 2:(p+1)){
    beta.hat[,j,] <- getBands(beta.samp[,j,], plot = plot, add = FALSE, x = tau.g, xlab = expression(tau), ylab = "Coefficient", bty = "n")
    if(plot) {
      title(main = plot.titles[j])
      abline(h = 0, lty = 2, col = 4)
    }
  }	
  # if(plot) suppressWarnings(par(cur.par,no.readonly = TRUE))  # return R parameters to pre-plot settings
  dimnames(beta.hat) <- list(tau=tau.g, beta=plot.titles, summary=c("b.lo", "b.med", "b.hi"))
  dimnames(beta.samp) <- list(tau=tau.g, beta=plot.titles, iter=1:length(ss))
  invisible(list(beta.samp = beta.samp, beta.est = beta.hat))
  
  mid.red <- which(tau.g == object$tau.g[mid])
  parametric.list <- rbind(beta.samp[mid.red, , ,drop=TRUE],
                           sigma = sigFn(pars[nknots * (p+1) + (p+1) + 1,ss], a.sig),	
                           nu = nuFn(pars[(nknots + 1) * (p+1) + 2,ss]))	
  dimnames(parametric.list)[[1]][1 + 0:p] <- c("Intercept", object$xnames)	
  gamsignu <- t(apply(parametric.list, 1, quantile, pr = c(0.5, 0.025, 0.975)))	
  dimnames(gamsignu)[[2]] <- c("Estimate", "Lo95%", "Up95%")	
  invisible(list(beta.samp = beta.samp, beta.est = beta.hat, parametric = gamsignu))
}


predict.qrjointtp <- function(object, x.test, location.test, tau.test, burn.perc = 0.5, nmc = 5e2){
  nsamp <- object$dim[11]; n <- object$dim[1]; m <- object$dim[5]
  nphi <- object$dim[8]; p <- object$dim[2]
  usamp <- matrix(object$usamp,ncol=n,byrow=T)
  parmgsamp <- matrix(object$parmgsamp, nrow = nsamp, byrow = TRUE)
  coeftp <- coef.qrjointtp(object, plot=FALSE, nmc = nmc, burn.perc = burn.perc, reduce = FALSE)
  ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 1)[-1]))
  
  n.test <- nrow(x.test)
  
  cte <- 2^(1-object$kappa)/gamma(object$kappa)
  
  tau.g <- object$tau.g
  len.reg.ix <- length(tau.g)
  tau.post <- array(NA, dim = c(nmc, n.test, len.reg.ix))
  tau.train.post <- array(NA, dim = c(nmc, n, len.reg.ix))
  
  
  cov.grid <- array(NA, dim = c(nphi, n, n))
  for(i in 1:nphi){
    cov.grid[i,,] <- covmat(phi = object$phi.grid[i], distance = object$distance, n = n, kernel = object$kernel, kappa = object$kappa, cte = cte)
  }
  
  distance.test <- t(apply(location.test, 1, function(u) apply(object$location, 1, function(v) sqrt(sum((u-v)^2)))))
  
  for(i in 1:nmc){
    psi <- plogis(object$lpsisamp[ss[i]])*18+2
    z.g <- qt(tau.g, df = psi)
    zsamp <- qt(usamp[ss[i],], df = psi)
    
    alpha <- exp(object$lkappasamp[ss[i]])/exp(parmgsamp[ss[i],(m+1)*(p+1)+1])
    sigma22 <- alpha*cov.grid[object$phisamp[ss[i]]+1,,]+(1-alpha)*diag(n)
    sigma22.inv <- chol2inv(chol(sigma22))
    dmaha <- t(zsamp)%*%sigma22.inv%*%zsamp
    
    for(j in 1:n.test){
      sigma11 <- 1
      sigma12 <- alpha*corfun(phi = object$phi.grid[object$phisamp[ss[i]]+1], distance = distance.test[j,], kernel = object$kernel, kappa = object$kappa, cte = cte)
      sigma12.22.inv <- sigma12%*%sigma22.inv
      cond.var <- as.numeric((psi + dmaha)*(sigma11 - tcrossprod(sigma12.22.inv, sigma12))/(psi + n))
      cond.mean <- c(sigma12.22.inv %*% zsamp)
      tau.post[i,j,] <- pt((z.g - cond.mean)/sqrt(cond.var), df = psi + n)
    }
  }
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  betahat <- function(k, n.loc, n.col, len.tau.comp){
    beta.hat <- array(NA, dim = c(n.loc, n.col, len.tau.comp))
    for(i in 1:n.loc){
      for(l in 1:n.col){
        beta.hat[i,l,] <- approx(x = tau.post[k,i,], y = coeftp$beta.samp[,l,k], xout = tau.test)$y
      }
    }
    return(beta.hat)
  }
  
  tau.len <- length(tau.test)
  beta.hat.list <- foreach(k=1:nmc) %dopar% {betahat(k, n.test, p+1, tau.len)}
  stopCluster(cl)
  
  beta.mean <- Reduce('+', beta.hat.list)/nmc
  
  cond.qt <- matrix(NA, nrow = n.test, ncol = tau.len)
  x.test.inter <- cbind(1, x.test)
  
  for(i in 1:n.test){
    for(k in 1:tau.len){
      cond.qt[i,k] <- sum(x.test.inter[i,]*beta.mean[i,,k])
    }
  }
  
  return(cond.qt)
}


#########################################################################


waictp <- function(object, burn.perc = 0.5, nmc = 200){
  nsamp <- object$dim[11]                      # number of iterations retains in sample
  pars <- matrix(object$parmgsamp, ncol = nsamp) # reorder parameters into npar x nsamp matrix
  ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 1)[-1])) # indices over which to summarizing (exclude burn; keep nmc samples)
  ss.len <- length(ss)
  
  dimpars <- object$dim; n <- dimpars[1]; p <- dimpars[2]; m <- dimpars[5]; ngrid <- dimpars[6]; nphi <- dimpars[8]; dimpars[9] <- ss.len; nparmg = (m+1)*(p+1) + 2;
  sm <- .C("DEV", pars = as.double(pars[,ss]), x = as.double(object$x), y = as.double(object$y), cens = as.integer(object$cens), wt = as.double(object$wt),
           shrink = as.integer(object$shrink), hyper = as.double(object$hyper), dim = as.integer(dimpars), gridmats = as.double(object$gridmats), tau.g = as.double(object$tau.g),
           devsamp = double(ss.len), llsamp = double(ss.len*n), pgsamp = double(ss.len*ngrid*(p+1)), rpsamp = double(ss.len*n), fbase.choice = as.integer(object$fbase.choice))
  deviance <- sm$devsamp
  ll <- matrix(sm$llsamp, ncol = ss.len)
  rp <- matrix(sm$rpsamp, ncol = ss.len)
  
  usamp <- matrix(object$usamp, ncol = nsamp)[,ss]
  
  phigrid <- matrix(object$phigrid, ncol = nphi)
  
  Lphigrid <- matrix(NA, nrow = nphi, ncol = n)
  Gphigrid <- Sigmagrid <- Sigmainvgrid <- array(NA, dim = c(nphi, n, n))
  
  for(i in 1:nphi){
    Lphigrid[i,] <- phigrid[1:n,i]
    Gphigrid[i,,] <- matrix(phigrid[-(1:n),i], ncol = n)
    Sigmagrid[i,,] <- Gphigrid[i,,]%*%diag(Lphigrid[i,])%*%t(Gphigrid[i,,])
    Sigmainvgrid[i,,] <- solve(Sigmagrid[i,,]) 
  }
  rhosamp <- exp(object$lkappasamp[ss])/exp(matrix(object$parmgsamp, ncol = nsamp)[nparmg - 1,ss])
  phisamp <- object$phisamp[ss] + 1
  psisamp <- plogis(object$lpsisamp[ss]) * 18 + 2
  wsamp <- epssamp <- matrix(NA, nrow = n, ncol = ss.len)
  
  
  for(i in 1:ss.len){
    Omega1 <- solve(rhosamp[i]*Sigmagrid[phisamp[i],,] + diag(1 - rhosamp[i], n))
    Omega1 <- (Omega1 + t(Omega1))/2
    zsamp <- qt(usamp[,i], df = psisamp[i])
    
    theta <- rgamma(1, shape = (psisamp[i] + n)/2, rate = (psisamp[i] + t(zsamp)%*%Omega1%*%zsamp)/2)
    
    Omega2 <- solve(Sigmainvgrid[phisamp[i],,] + diag(rhosamp[i]/(1 - rhosamp[i]), n))
    Omega2 <- (Omega2 + t(Omega2))/2
    
    wsamp[,i] <- rmvnorm(1, mean = sqrt(rhosamp[i]*theta)/(1 - rhosamp[i])*Omega2%*%zsamp, sigma = Omega2)
    epssamp[,i] <- (zsamp - sqrt(rhosamp[i]/theta)*wsamp[,i])/sqrt((1 - rhosamp[i])/theta)
    ll[,i] <- ll[,i] - log(sqrt((1 - rhosamp[i])/theta)*dt(zsamp, df = psisamp[i])/dnorm(epssamp[,i]))
  }
  
  fit.waic <- waic(ll)
  return(fit.waic)
}


# Function to construct beta covariate coefficient functionals from posterior outputs

estFn <- function(par, x, y, gridmats, L, mid, nknots, ngrid, a.kap, a.sig, tau.g, reg.ix, reduce = TRUE, x.ce = 0, x.sc = 1, base.bundle){
  
  n <- length(y); p <- ncol(x)
  wKnot <- matrix(par[1:(nknots*(p+1))], nrow = nknots)
  w0PP <- ppFn0(wKnot[,1], gridmats, L, nknots, ngrid)
  w0 <- w0PP$w
  wPP <- apply(wKnot[,-1,drop=FALSE], 2, ppFn, gridmats = gridmats, L = L, nknots = nknots, ngrid = ngrid, a.kap = a.kap)
  wMat <- matrix(sapply(wPP, extract, vn = "w"), ncol = p)
  
  zeta0.dot <- exp(shrinkFn(p) * (w0 - max(w0)))          # subract max for numerical stability
  zeta0 <- trape(zeta0.dot[-c(1,L)], tau.g[-c(1,L)], L-2) # take integral of derivative to get original function.
  zeta0.tot <- zeta0[L-2] 
  zeta0 <- c(0, tau.g[2] + (tau.g[L-1]-tau.g[2])*zeta0 / zeta0.tot, 1)
  zeta0.dot <- (tau.g[L-1]-tau.g[2])*zeta0.dot / zeta0.tot
  zeta0.dot[c(1,L)] <- 0
  zeta0.ticks <- pmin(L-1, pmax(1, sapply(zeta0, function(u) sum(tau.g <= u))))
  zeta0.dists <- (zeta0 - tau.g[zeta0.ticks]) / (tau.g[zeta0.ticks+1] - tau.g[zeta0.ticks])
  vMat <- apply(wMat, 2, transform.grid, ticks = zeta0.ticks, dists = zeta0.dists)
  
  reach <- nknots*(p+1)
  gam0 <- par[reach + 1]; reach <- reach + 1
  gam <- par[reach + 1:p]; reach <- reach + p
  sigma <- sigFn(par[reach + 1], a.sig); reach <- reach + 1
  nu <- nuFn(par[reach + 1]);
  
  # Intercept quantiles
  b0dot <- sigma * base.bundle$q0(zeta0, nu) * zeta0.dot
  beta0.hat <- rep(NA, L)
  beta0.hat[mid:L] <- gam0 + trape(b0dot[mid:L], tau.g[mid:L], L - mid + 1)
  beta0.hat[mid:1] <- gam0 + trape(b0dot[mid:1], tau.g[mid:1], mid)
  
  # Working towards other covariate coefficient quantiles
  # First four lines get 'aX' shrinkage factor
  vNorm <- sqrt(rowSums(vMat^2))
  a <- tcrossprod(vMat, x)
  aX <- apply(-a, 1, max)/vNorm
  aX[is.nan(aX)] <- Inf
  aTilde <- vMat / (aX * sqrt(1 + vNorm^2))
  ab0 <- b0dot * aTilde
  
  beta.hat <- kronecker(rep(1,L), t(gam))
  beta.hat[mid:L,] <- beta.hat[mid:L,] + apply(ab0[mid:L,,drop=FALSE], 2, trape, h = tau.g[mid:L], len = L - mid + 1)
  beta.hat[mid:1,] <- beta.hat[mid:1,] + apply(ab0[mid:1,,drop=FALSE], 2, trape, h = tau.g[mid:1], len = mid)
  beta.hat <- beta.hat / x.sc
  beta0.hat <- beta0.hat - rowSums(beta.hat * x.ce)
  betas <- cbind(beta0.hat, beta.hat)
  if(reduce) betas <- betas[reg.ix,,drop = FALSE]
  return(betas)
}


chull.center <- function (x, maxEPts = ncol(x) + 1, plot = FALSE){
  sx <- as.matrix(apply(x, 2, function(s) punif(s, min(s), max(s))))
  dd <- rowSums(scale(sx)^2)
  ix.dd <- order(dd, decreasing = TRUE)
  sx <- sx[ix.dd, , drop = FALSE]
  x.chol <- inchol(sx, maxiter = maxEPts)   # uses default rbfdot (radial basis kernal function, "Gaussian")
  ix.epts <- ix.dd[pivots(x.chol)]          # picks off extreme points
  x.choose <- x[ix.epts, , drop = FALSE]
  xCent <- as.numeric(colMeans(x.choose))   # takes means of extreme points
  attr(xCent, "EPts") <- ix.epts            # returns list of extrememum points
  if (plot) {
    n <- nrow(x)
    p <- ncol(x)
    xjit <- x + matrix(rnorm(n * p), n, p) %*% diag(0.05 * apply(x, 2, sd), p)
    xmean <- colMeans(x)
    x.ept <- x[ix.epts, ]
    M <- choose(p, 2)
    xnames <- dimnames(x)[[2]]
    if (is.null(xnames)) xnames <- paste("x", 1:p, sep = "")
    par(mfrow = c(ceiling(M/ceiling(sqrt(M))), ceiling(sqrt(M))), mar = c(2, 3, 0, 0) + 0.1)
    xmax <- apply(x, 2, max)
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        plot(x[, i], x[, j], col = "gray", cex = 1, ann = FALSE, ty = "n", axes = FALSE, bty = "l")
        title(xlab = xnames[i], line = 0.3)
        title(ylab = xnames[j], line = 0.3)
        ept <- chull(x[, i], x[, j])
        polygon(x[ept, i], x[ept, j], col = gray(0.9), border = "white")
        points(xjit[, i], xjit[, j], pch = ".", col = gray(0.6))
        points(xmean[i], xmean[j], col = gray(0), pch = 17, cex = 1)
        points(xCent[i], xCent[j], col = gray(0), pch = 1, cex = 2)
        points(x.ept[, i], x.ept[, j], col = gray(.3), pch = 10, cex = 1.5)
      }
    }
  }
  return(xCent)
}


waic <- function(logliks, print = TRUE){
  lppd <- sum(apply(logliks, 1, logmean))
  p.waic.1 <- 2 * lppd - 2 * sum(apply(logliks, 1, mean))
  p.waic.2 <- sum(apply(logliks, 1, var))
  waic.1 <- -2 * lppd + 2 * p.waic.1
  waic.2 <- -2 * lppd + 2 * p.waic.2
  if(print) cat("WAIC.1 =", round(waic.1, 2), ", WAIC.2 =", round(waic.2, 2), "\n")
  invisible(c(WAIC1 = waic.1, WAIC2 = waic.2))
}

# Returns W-tilde_0 (not explictly mentioned in the paper)
# Using t-distribution with three degrees of freedom as prior distribution
# get L-dimensional vector that is needed for likelihood evaluation.

ppFn0 <- function(w.knot, gridmats, L, nknots, ngrid){
  w.grid <- matrix(NA, L, ngrid)
  lpost.grid <- rep(NA, ngrid)
  for(i in 1:ngrid){
    A <- matrix(gridmats[1:(L*nknots),i], nrow = nknots)
    R <- matrix(gridmats[L*nknots + 1:(nknots*nknots),i], nrow = nknots)
    r <- sum.sq(backsolve(R, w.knot, transpose = TRUE))
    w.grid[,i] <- colSums(A * w.knot)
    lpost.grid[i] <- -(0.5*nknots+0.1)*log1p(0.5*r/0.1) - gridmats[nknots*(L+nknots)+1,i] + gridmats[nknots*(L+nknots)+2,i]		
  }
  lpost.sum <- logsum(lpost.grid)
  post.grid <- exp(lpost.grid - lpost.sum)
  w <- c(w.grid %*% post.grid)
  return(list(w = w, lpost.sum = lpost.sum))
}

# Returns W-tilde_j used in the paper, L-dimensional vector needed for likelihood eval.
# Also returns the log posterior sum.

ppFn <- function(w.knot, gridmats, L, nknots, ngrid, a.kap){
  w.grid <- matrix(NA, L, ngrid)
  lpost.grid <- rep(NA, ngrid)
  for(i in 1:ngrid){
    A <- matrix(gridmats[1:(L*nknots),i], nrow = nknots)
    R <- matrix(gridmats[L*nknots + 1:(nknots*nknots),i], nrow = nknots)
    r <- sum.sq(backsolve(R, w.knot, transpose = TRUE))
    w.grid[,i] <- colSums(A * w.knot)
    lpost.grid[i] <- (logsum(-(nknots/2+a.kap[1,])*log1p(0.5*r/ a.kap[2,]) + a.kap[3,] + lgamma(a.kap[1,]+nknots/2)-lgamma(a.kap[1,])-.5*nknots*log(a.kap[2,]))
                      - gridmats[nknots*(L+nknots)+1,i] + gridmats[nknots*(L+nknots)+2,i])		
  }
  lpost.sum <- logsum(lpost.grid)
  post.grid <- exp(lpost.grid - lpost.sum)
  w <- c(w.grid %*% post.grid)
  return(list(w = w, lpost.sum = lpost.sum))
}


lamFn <- function(prox) return(sqrt(-100*log(prox)))
nuFn <- function(z) return(0.5 + 5.5*exp(z/2)) 
nuFn.inv <- function(nu) return(2*log((nu - 0.5)/5.5))
sigFn <- function(z, a.sig) return(exp(z/2)) 
sigFn.inv <- function(s, a.sig) return(2 * log(s))
unitFn <- function(u) return(pmin(1 - 1e-10, pmax(1e-10, u)))

sum.sq <- function(x) return(sum(x^2))
extract <- function(lo, vn) return(lo[[vn]])       # pull out the "vn"th list item
logmean <- function(lx) return(max(lx) + log(mean(exp(lx - max(lx)))))
logsum <- function(lx) return(logmean(lx) + log(length(lx)))
shrinkFn <- function(x) return(1) ##(1/(1 + log(x)))
trape <- function(x, h, len = length(x)) return(c(0, cumsum(.5 * (x[-1] + x[-len]) * (h[-1] - h[-len]))))

getBands <- function(b, col = 2, lwd = 1, plot = TRUE, add = FALSE, x = seq(0,1,len=nrow(b)), remove.edges = TRUE, ...){
  colRGB <- col2rgb(col)/255
  colTrans <- rgb(colRGB[1], colRGB[2], colRGB[3], alpha = 0.2)
  b.med <- apply(b, 1, quantile, pr = .5)
  b.lo <- apply(b, 1, quantile, pr = .025)
  b.hi <- apply(b, 1, quantile, pr = 1 - .025)
  L <- nrow(b)
  ss <- 1:L; ss.rev <- L:1
  if(remove.edges){
    ss <- 2:(L-1); ss.rev <- (L-1):2
  }
  if(plot){
    if(!add) plot(x[ss], b.med[ss], ty = "n", ylim = range(c(b.lo[ss], b.hi[ss])), ...)  
    polygon(x[c(ss, ss.rev)], c(b.lo[ss], b.hi[ss.rev]), col = colTrans, border = colTrans)
    lines(x[ss], b.med[ss], col = col, lwd = lwd)  
  }
  invisible(cbind(b.lo, b.med, b.hi))
}

# Function to calculate kulback-liebler divergence between two GPS with alternate
# lambdas, detailed in appendix B.2
# Uses formula for kullback liebler divergence of two multiva gauassian dist'ns,
# each of dimension k-knots
# solve(K2,K1) does K2^{-1}%*%K1... and sum/diag around it gives determinant

klGP <- function(lam1, lam2, nknots = 11){
  tau <- seq(0, 1, len = nknots)
  dd <- outer(tau, tau, "-")^2
  K1 <- exp(-lam1^2 * dd); diag(K1) <- 1 + 1e-10; R1 <- chol(K1); log.detR1 <- sum(log(diag(R1)))
  K2 <- exp(-lam2^2 * dd); diag(K2) <- 1 + 1e-10; R2 <- chol(K2); log.detR2 <- sum(log(diag(R2)))
  return(log.detR2-log.detR1 - 0.5 * (nknots - sum(diag(solve(K2, K1)))))
}

# Function to create non-evenly-spaced grid for lambda discretized grid points
# Involves using kulback-liebler divergence.  Makes sure that prior remains sufficiently
# overlapped for neighboring lambda values, preventing (hopefully) poor mixing of
# Markov chain sampler.

proxFn <- function(prox.Max, prox.Min, kl.step = 1){
  prox.grid <- prox.Max
  j <- 1
  while(prox.grid[j] > prox.Min){
    prox1 <- prox.grid[j]
    prox2 <- prox.Min
    kk <- klGP(lamFn(prox1), lamFn(prox2))
    while(kk > kl.step){
      prox2 <- (prox1 + prox2)/2
      kk <- klGP(lamFn(prox1), lamFn(prox2))
    }
    j <- j + 1
    prox.grid <- c(prox.grid, prox2)
  }
  return(prox.grid)
}

transform.grid <- function(w, ticks, dists){
  return((1-dists) * w[ticks] + dists * w[ticks+1])
}

