library(plyr)

qrjointtp <- function(x, y, location = NULL, distance = NULL, nphi, kappa, kernel = c("se", "matern"), nsamp = 1e3, thin = 10, cens = NULL,
                       wt = NULL, incr = 0.01, par = "prior", nknots = 6,
                       hyper = list(sig = c(.1,.1), lam = c(6,4), kap = c(0.1,0.1,1)),
                       shrink = FALSE, prox.range = c(.2,.95), acpt.target = 0.15, acpt.hptarget = 0.15,
                       acpt.lpsitarget = 0.15,
                       ref.size = 3, ref.hpsize = 3, ref.lpsisize = 3, blocking = "std5", temp = 1, expo = 2,
                       blocks.mu, blocks.S, fix.nu = FALSE, fbase = c("t", "logistic", "unif"), verbose = TRUE){
  
  # Set up base functions
  fbase.choice <- match(fbase[1], c("t", "logistic", "unif"))
  if(is.na(fbase.choice)) stop("Only 't', 'logistic' or 'unif' is allowed for the choice of fbase")
  if(fbase.choice == 1){
    q0 <- function(u, nu = Inf) return(1 / (dt(qt(unitFn(u), df = nu), df = nu) * qt(.9, df = nu)))
    Q0 <- function(u, nu = Inf) return(qt(unitFn(u), df = nu) / qt(.9, df = nu))
    F0 <- function(x, nu = Inf) return(pt(x*qt(.9, df = nu), df = nu))
  } else if(fbase.choice == 2){
    fix.nu <- 1
    q0 <- function(u, nu = Inf) return(1 / dlogis(qlogis(unitFn(u))))
    Q0 <- function(u, nu = Inf) return(qlogis(unitFn(u)))
    F0 <- function(x, nu = Inf) return(plogis(x))
  } else {
    fix.nu <- 1
    q0 <- function(u, nu = Inf) return(1 / (dunif(qunif(u, -1,1), -1,1)))
    Q0 <- function(u, nu = Inf) return(qunif(u, -1,1))
    F0 <- function(x, nu = Inf) return(punif(x, -1,1))
  }
  base.bundle <- list(q0 = q0, Q0 = Q0, F0 = F0)
  
  n <- nrow(x)
  p <- ncol(x)
  x <- scale(x, chull.center(x))
  
  if(is.null(cens)) cens <- rep(0,n) else{
    cens <- as.vector(cens)
    if(!is.numeric(cens))
      stop("'cens' must be a numeric vector")
    
    if (length(cens) != n) 
      stop(gettextf("number of cens is %d, should equal %d (number of observations)", 
                    length(cens), n), domain = NA)}
  
  if(is.null(wt)) wt <- rep(1,n) else{
    wt <- as.vector(wt)
    if(!is.numeric(wt))
      stop("'wt' must be a numeric vector")
    
    if (length(wt) != n) 
      stop(gettextf("number of wt is %d, should equal %d (number of observations)", 
                    length(wt), n), domain = NA)}
  
  # tau.g:  sequence of tau probabilities to be spead 'incr' apart, 
  #         supplemented at upper and lower regions by additionally fine grid of taus
  # L:      number of tau locations
  # mid:    indexes prob .5, (e.g. median) location in tau grid
  # reg.ix: indices of tau belonging to non-supplemented grid	
  
  Ltail <- ceiling(2*log(n,2) + log(incr,2))
  if(Ltail > 0){
    tau.tail <- incr / 2^(Ltail:1)
    tau.g <- c(0, tau.tail, seq(incr, 1 - incr, incr), 1 - tau.tail[Ltail:1], 1)
    L <- length(tau.g); mid <- which.min(abs(tau.g - 0.5)); reg.ix <- (1:L)[-c(1 + 1:Ltail, L - 1:Ltail)]
  } else {
    tau.g <- seq(0, 1, incr)
    L <- length(tau.g); mid <- which.min(abs(tau.g - 0.5)); reg.ix <- (1:L)
  }
  
  tau.kb <- seq(0,1,len = nknots)
  tau.k <- tau.kb
  
  # set priors to defaults if not specified
  a.sig <- hyper$sig; if(is.null(a.sig)) a.sig <- c(.1, .1)
  a.lam <- hyper$lam; if(is.null(a.lam)) a.lam <- c(6, 4)
  a.kap <- hyper$kap; if(is.null(a.kap)) a.kap <- c(0.1, 0.1, 1); a.kap <- matrix(a.kap, nrow = 3); nkap <- ncol(a.kap); a.kap[3,] <- log(a.kap[3,])
  hyper.reduced <- c(a.sig, c(a.kap))
  
  # Create grid-discretized lambdas; retain sufficient overlap to get good mixing.
  # User specifies 1) range for reasonable correlations and 2) parameters for 
  # beta distribution dictating the correlation probabilities.
  # Code translates to lambda scale & gives discretized probs
  prox.grid <- proxFn(max(prox.range), min(prox.range), 0.5)
  ngrid <- length(prox.grid)
  lamsq.grid <- lamFn(prox.grid)^2
  prior.grid <- -diff(pbeta(c(1, (prox.grid[-1] + prox.grid[-ngrid])/2, 0), a.lam[1], a.lam[2]))
  lp.grid <- log(prior.grid)  # log probabilities of the lambda prior
  
  d.kg <- abs(outer(tau.k, tau.g, "-"))^expo	# dist between orig tau grid and reduced tau grid
  d.kk <- abs(outer(tau.k, tau.k, "-"))^expo  # dist between points on reduced tau grid
  gridmats <- matrix(NA, nknots*(L + nknots)+2, ngrid) # setup to hold stuff; (long) x dim of lambda grid
  K0 <- 0
  t1 <- Sys.time()
  for(i in 1:ngrid){
    K.grid <- exp(-lamsq.grid[i] * d.kg); K.knot <- exp(-lamsq.grid[i] * d.kk);	diag(K.knot) <- 1 + 1e-10	
    R.knot <- chol(K.knot); A.knot <- solve(K.knot, K.grid)   # Chol decomp of covariance (not yet inverted)
    # gridmats: Concatenation of first A_0g (L x m), then R_0g (m x m), 
    # then log det of R, then log lambda prior.  Done for each lambda separately	
    gridmats[,i] <- c(c(A.knot), c(R.knot), sum(log(diag(R.knot))), lp.grid[i])	
    # K0: Prior correlation, weighted by prior on different lambda values
    K0 <- K0 + prior.grid[i] * K.knot
  }
  t2 <- Sys.time()
  #cat("Matrix calculation time per 1e3 iterations =", round(1e3 * as.numeric(t2 - t1), 2), "\n")
  
  ## create list of parameters to be passed to c
  # n: number of observations
  # p: number of covariates in design matrix, not inluding intercept
  # L: number of tau grid points
  # mid: index of tau in grid acting as median.  (indexed as one fewer in C than R)
  # nknots: number of knots in interpolation approximation of each w_j
  # ngrid: number of points in lambda grid/discretized prior (lambda dictates
  #        length-scale in GP, ie how related different taus are to eachother)
  
  niter <- nsamp * thin
  dimpars <- c(n, p, L, mid - 1, nknots, ngrid, ncol(a.kap), nphi, niter, thin, nsamp)
  nparmg <- (nknots+1) * (p+1) + 2
  # First supported option: using prior to initialize MC iteration
  if(par[1] == "prior") {
    parmg <- rep(0, nparmg)   # vector to hold all parameters
    if(fix.nu) parmg[nparmg] <- nuFn.inv(fix.nu) # nu in last slot of par
    
    # Get rq beta coefficients for each tau in grid. Then use 5th degree b-spline
    # to estimate curve for each beta coefficient
    # "dither" slightly perturbs responses in order to avoid potential degeneracy
    # of dual simplex algorithm (present with large number of ties in y) and "hanging"
    # of the rq Fortran code
    beta.rq <- sapply(tau.g, function(a) return(coef(suppressWarnings(rq(dither(y) ~ x, tau = a, weights = wt)))))   # THIS COULD BE DONE WITH tau=tau.g
    v <- bs(tau.g, df = 5)
    
    # over tau and per coefficient get smoothed fits through data (ie independently estimated quantiles);
    # ultimate goal: reasonable estimates at median, at delta and at 1-delta
    rq.lm <- apply(beta.rq, 1, function(z) return(coef(lm(z ~ v))))
    
    # Combine b-spline estimates to get estimate of coefficient when tau
    # is delta (lower quantile) and when tau is 1 - delta (upper quantile)
    delta <- tau.g[2]
    tau.0 <- tau.g[mid]
    rq.tau0 <- c(c(1, predict(v, tau.0)) %*% rq.lm)
    rq.delta <- c(c(1, predict(v, delta)) %*% rq.lm)
    rq.deltac <- c(c(1, predict(v, 1 - delta)) %*% rq.lm)
    
    parmg[nknots*(p+1) + 1:(p+1)] <- as.numeric(rq.tau0) # median coefficient quantiles  placed in p + 1 positions just before last 2
    sigma <- 1
    parmg[(nknots+1)*(p+1) + 1] <- sigFn.inv(sigma, a.sig) # prior for sigma.  
    
    for(i in 1:(p+1)){
      kapsq <- sum(exp(a.kap[3,]) * (a.kap[2,] / a.kap[1,]))
      lam.ix <- sample(length(lamsq.grid), 1, prob = prior.grid)
      R <- matrix(gridmats[L*nknots + 1:(nknots*nknots),lam.ix], nknots, nknots)
      z <- sqrt(kapsq) * c(crossprod(R, rnorm(nknots)))
      parmg[(i - 1) * nknots + 1:nknots] <- z - mean(z) # centered and placed in sets of length nknots one after another for p+1 covariates
    }
    
    # get betas associated with these prior params & corresponding quantile estimate (centered around median)
    beta.hat <- estFn(parmg, x, y, gridmats, L, mid, nknots, ngrid, a.kap, a.sig, tau.g, reg.ix, FALSE, base.bundle = base.bundle)
    qhat <- tcrossprod(cbind(1, x), beta.hat)
    
    infl <- max(max((y - qhat[,mid])/(qhat[,ncol(qhat) - 1] - qhat[,mid])), max((qhat[,mid] - y)/(qhat[,mid] - qhat[,2])))
    oo <- .C("INIT", par = as.double(parmg), x = as.double(x), y = as.double(y), cens = as.integer(cens), wt = as.double(wt),
             shrink = as.integer(shrink), hyper = as.double(hyper.reduced), dim = as.integer(dimpars), gridpars = as.double(gridmats),
             tau.g = as.double(tau.g), siglim = as.double(sigFn.inv(c(1.0 * infl * sigma, 10 * infl * sigma), a.sig)),
             fbase.choice = as.integer(fbase.choice))
    
    parmg <- oo$par
  } else if (par[1] == "RQ"){
    
    # Initialize at a model space approximation of the estimates from rq
    parmg <- rep(0, (nknots+1) * (p+1) + 2)
    delta <- tau.g[2]     # smallest but one
    tau.0 <- tau.g[mid]   # median or mid-most value of taus evaluated
    
    epsilon <- 0.1 * min(diff(sort(tau.k)))
    tau.knot.plus <- pmin(tau.k + epsilon, 1)
    tau.knot.minus <- pmax(tau.k - epsilon, 0)
    
    beta.rq <- sapply(c(tau.0, delta, 1 - delta, tau.knot.plus, tau.knot.minus), function(a) return(coef(suppressWarnings(rq(y~ x, tau = a, weights = wt)))))
    
    rq.tau0 <- c(beta.rq[,1])       # estimate coef median MLEs 
    rq.delta <- c(beta.rq[,2])     # estimate coef upper quantile MLEs
    rq.deltac <- c(beta.rq[,3])  # estimate coef lower quantile MLEs
    
    parmg[nknots*(p+1) + 1:(p+1)] <- as.numeric(rq.tau0)
    nu <- ifelse(fix.nu, fix.nu, nuFn(0))
    # Solving for sigma in zeta(tau) - F_0( (beta0(delta) - gam0)/sigma)
    sigma <- max((rq.delta[1] - rq.tau0[1]) / (Q0(delta, nu) - rq.tau0[1]), (rq.deltac[1] - rq.tau0[1]) / (Q0(1 - delta, nu) - rq.tau0[1]))
    # sigma <- min((rq.delta[1] - rq.tau0[1]) / Q0(delta, nu), (rq.deltac[1] - rq.tau0[1]) / Q0(1 - delta, nu))
    parmg[(nknots+1)*(p+1) + 1]  <- sigFn.inv(sigma, a.sig)
    
    beta.rq.plus <- t(beta.rq[,3+1:nknots])
    beta.rq.minus <- t(beta.rq[,3+nknots+1:nknots])
    
    
    
    zeta0.plus <- F0((beta.rq.plus[,1] - rq.tau0[1]) / sigma + rq.tau0[1], nu)
    zeta0.minus <- F0((beta.rq.minus[,1] - rq.tau0[1]) / sigma + rq.tau0[1], nu)
    zeta0.dot.knot <- (zeta0.plus - zeta0.minus) / (tau.knot.plus - tau.knot.minus)
    w0.knot <- log(pmax(epsilon, zeta0.dot.knot)) / shrinkFn(p)
    w0.knot <- (w0.knot - mean(w0.knot))
    
    w0PP <- ppFn0(w0.knot, gridmats, L, nknots, ngrid)
    w0 <- w0PP$w
    
    zeta0.dot <- exp(shrinkFn(p) * (w0 - max(w0)))
    zeta0 <- trape(zeta0.dot[-c(1,L)], tau.g[-c(1,L)], L-2)
    zeta0.tot <- zeta0[L-2]
    zeta0 <- c(0, tau.g[2] + (tau.g[L-1]-tau.g[2])*zeta0 / zeta0.tot, 1)
    zeta0.dot <- (tau.g[L-1]-tau.g[2])*zeta0.dot / zeta0.tot
    zeta0.dot[c(1,L)] <- 0
    
    parmg[1:nknots] <- w0.knot
    beta0.dot <- sigma * q0(zeta0, nu) * zeta0.dot
    
    tilt.knot <- tau.g[tilt.ix <- sapply(tau.k, function(a) which.min(abs(a - zeta0)))]
    
    tau.knot.plus <- pmin(tilt.knot + epsilon, 1)
    tau.knot.minus <- pmax(tilt.knot - epsilon, 0)
    beta.rq.plus <- t(sapply(tau.knot.plus, function(a) return(coef(suppressWarnings(rq(y~ x, tau = a, weights = wt))))))
    beta.rq.minus <- t(sapply(tau.knot.minus, function(a) return(coef(suppressWarnings(rq(y~ x, tau = a, weights = wt))))))
    
    beta.dot.knot <- (beta.rq.plus[,-1,drop=FALSE] - beta.rq.minus[,-1,drop=FALSE]) /  (tau.knot.plus - tau.knot.minus)
    
    parmg[nknots + 1:(nknots*p)] <- c(beta.dot.knot)
    beta.hat <- estFn(parmg, x, y, gridmats, L, mid, nknots, ngrid, a.kap, a.sig, tau.g, reg.ix, FALSE, base.bundle = base.bundle)
    qhat <- tcrossprod(cbind(1, x), beta.hat)
    infl <- max(max((y - qhat[,mid])/(qhat[,ncol(qhat) - 1] - qhat[,mid])), max((qhat[,mid] - y)/(qhat[,mid] - qhat[,2])))
    parmg[(nknots+1)*(p+1) + 1]  <- sigFn.inv(sigma, a.sig)
    oo <- .C("INIT", par = as.double(parmg), x = as.double(x), y = as.double(y), cens = as.integer(cens),
             wt = as.double(wt), shrink = as.integer(shrink), hyper = as.double(hyper.reduced), dim = as.integer(dimpars), gridpars = as.double(gridmats),
             tau.g = as.double(tau.g), siglim = as.double(sigFn.inv(c(1.0 * infl * sigma, 10.0 * infl * sigma), a.sig)), fbase.choice = as.integer(fbase.choice))

    parmg <- oo$par
  }
  
  # Choose a blocking scheme for updates.
  if(blocking == "single"){                                                # ONE BLOCK
    blocks <- list(rep(TRUE, nparmg))                                        # Update all simultaneously
  } else if(blocking == "single2"){                                        # TWO BLOCKS
    blocks <- list(rep(TRUE, nparmg), rep(FALSE, nparmg))                      # First, update all simultaneously
    blocks[[2]][nknots*(p+1) + 1:(p+3)] <- TRUE                            # Second, update lambda_nots, sigma, and nu simultaneously
  } else if(blocking == "single3"){                                        # THREE BLOCKS
    blocks <- list(rep(TRUE, nparmg), rep(FALSE, nparmg), rep(FALSE, nparmg))    # First, update all simultaneously
    blocks[[2]][nknots*(p+1) + 1:(p+1)] <- TRUE                            # Second, update all lambda_nots
    blocks[[3]][(nknots+1)*(p+1) + 1:2] <- TRUE                            # Last, update sigma and nu
  } else if(blocking == "std0"){                                           # P + 1 BLOCKS
    blocks <- replicate(p + 1, rep(FALSE, nparmg), simplify = FALSE)         # Update everything related to covariate (GP, lamda_not) simulatneous with sigma, nu
    for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots, nknots*(p+1) + i + 1, (nknots+1)*(p+1) + 1:2)] <- TRUE
  } else if(blocking == "std1"){                                           # P + 2 BLOCKS
    blocks <- replicate(p + 2, rep(FALSE, nparmg), simplify = FALSE)         # First, Lowrank GPS + sigma, nu for each covariate
    for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots, nknots*(p+1) + i + 1, (nknots+1)*(p+1) + 1:2)] <- TRUE
    blocks[[p + 2]][nknots*(p+1) + 1:(p+3)] <- TRUE                        # Then all covariate medians, sigma, nu simultaneously 
  } else if(blocking == "std2"){                                           # P + 3 BLOCKS
    blocks <- replicate(p + 3, rep(FALSE, nparmg), simplify = FALSE)         # Just like P+2, with additional update that is only simga and nu
    for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots, nknots*(p+1) + i + 1, (nknots+1)*(p+1) + 1:2)] <- TRUE
    blocks[[p + 2]][nknots*(p+1) + 1:(p+1)] <- TRUE
    blocks[[p + 3]][(nknots+1)*(p+1) + 1:2] <- TRUE
  } else if(blocking == "std3"){                                           # ALSO P + 3 BLOCKS
    blocks <- replicate(p + 3, rep(FALSE, nparmg), simplify = FALSE)
    for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots)] <- TRUE        # Update GPs related to covariates
    blocks[[p + 2]][nknots*(p+1) + 1:(p+1)] <- TRUE                        # Update covariate medians
    blocks[[p + 3]][(nknots+1)*(p+1) + 1:2] <- TRUE                        # Lastly, update sigma and nu
  } else if(blocking == "std4"){
    blocks <- replicate(p + 3, rep(FALSE, nparmg), simplify = FALSE)         # ALSO P + 3 BLOCKS
    for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots, nknots*(p+1) + i + 1)] <- TRUE  # Update GPs and medians
    blocks[[p + 2]][nknots*(p+1) + 1:(p+1)] <- TRUE                        # Update medians
    blocks[[p + 3]][(nknots+1)*(p+1) + 1:2] <- TRUE                        # Update sigma and nu
  } else if(blocking == "std5"){
    blocks <- replicate(p + 4, rep(FALSE, nparmg), simplify = FALSE)        # Same as previous with additional update to ALL parameters
    for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots, nknots*(p+1) + i + 1)] <- TRUE
    blocks[[p + 2]][nknots*(p+1) + 1:(p+1)] <- TRUE
    blocks[[p + 3]][(nknots+1)*(p+1) + 1:2] <- TRUE
    blocks[[p + 4]][1:nparmg] <- TRUE
  } else {                                                                # Univariate updates
    blocks <- replicate(nparmg, rep(FALSE, nparmg), simplify = FALSE)
    for(i in 1:nparmg) blocks[[i]][i] <- TRUE
  }
  
  nblocks <- length(blocks)
  if(fix.nu) for(j in 1:nblocks) blocks[[j]][(nknots+1) * (p+1) + 2] <- FALSE
  
  blocks.ix <- c(unlist(lapply(blocks, which))) - 1                       # Index of locations of updates
  blocks.size <- sapply(blocks, sum)  # Size of MVN in each block update
  if(missing(blocks.mu)) blocks.mu <- rep(0, sum(blocks.size))
  if(missing(blocks.S)){
    sig.nu <- c(TRUE, !as.logical(fix.nu))
    blocks.S <- lapply(blocks.size, function(q) diag(1, q))
    if(substr(blocking, 1, 3) == "std"){
      for(i in 1:(p+1)) blocks.S[[i]][1:nknots, 1:nknots] <- K0
      if(as.numeric(substr(blocking, 4,5)) > 1){
        blocks.S[[p + 2]] <- summary(suppressWarnings(rq(dither(y) ~ x, tau = 0.5, weights = wt)), se = "boot", cov = TRUE)$cov
        blocks.S[[p + 3]] <- matrix(c(1, 0, 0, .1), 2, 2)[sig.nu,sig.nu]
      }
      if(as.numeric(substr(blocking, 4,5)) == 5){
        slist <- list(); length(slist) <- p + 3
        for(i in 1:(p+1)) slist[[i]] <- K0
        slist[[p+2]] <- summary(suppressWarnings(rq(dither(y) ~ x, tau = 0.5, weights = wt)), se = "boot", cov = TRUE)$cov
        slist[[p+3]] <- matrix(c(1, 0, 0, .1), 2, 2)[sig.nu,sig.nu]
        blocks.S[[p + 4]] <- as.matrix(bdiag(slist))
      }
    }
    blocks.S <- unlist(blocks.S)
  }
  
  imcmc.par <- c(nblocks, ref.size, verbose, max(10, niter/1e4), ref.hpsize, 0, ref.lpsisize, 0, rep(0, nblocks))
  dmcmc.par <- c(temp, 0.999, acpt.hptarget, 2.38, acpt.lpsitarget, 2.38, rep(acpt.target, nblocks), 2.38 / sqrt(blocks.size))
  
  # Updated in C code
  # parmgsamp:     holds all marginal parameters from all iterations in one long vector
  # acptsamp:      holds all acceptance rates for each block at each iteration in one long vector
  # lpsamp:        holds all evaluations of log posterior at each iteration 
  # lkappasamp:    holds log copula parameter sigma_s^2 at each iteration 
  # lpsisamp:      holds log psi at each iteration in one long vector
  # phisamp:       holds copula parameter phi at each iteration
  # usamp:         holds all latent quantile levels at each iteration
  
  
  kernel.choice <- match(kernel[1], c("se", "matern"))
  
  cte <- 2^(1-kappa)/gamma(kappa)

  phi.lb <- uniroot(function(u) sapply(u, kernel[1], x = max(distance)/4, kappa = kappa, cte = cte) - 0.05, interval = c(1e-4, 10*max(distance)), tol = 1e-10)$root
  phi.ub <- 3*phi.lb
  phi.grid <- seq(phi.lb, phi.ub, length.out = nphi)
  cov.grid <- lapply(phi.grid, covmat, n = n, distance = distance, kernel = kernel.choice, kappa = kappa, cte = cte)
  eigen.grid <- lapply(cov.grid, eigen)
  
  rho <- 0.5
  phi <- floor(nphi/2)
  psi <- 10

  theta <- iota <- 0
  eta <- delta <- 1
  
  
  
  tm.c <- system.time(oo <- .C("BJQRtp", par = as.double(parmg), phi = as.integer(phi), rho = as.double(rho), psi = as.double(psi), 
                               phigrid = as.double(sapply(eigen.grid, function(u) cbind(u$val, u$vec))), 
                               x = as.double(x), y = as.double(y), cens = as.integer(cens), wt = as.double(wt),
                               shrink = as.integer(shrink), hyper = as.double(hyper.reduced), dim = as.integer(dimpars), 
                               gridmats = as.double(gridmats), tau.g = as.double(tau.g), muV = as.double(blocks.mu), SV = as.double(blocks.S), 
                               thetaVar = as.double(theta), etaVar = as.double(eta), iotaVar = as.double(iota), deltaVar = as.double(delta),
                               blocks = as.integer(blocks.ix), blocks.size = as.integer(blocks.size), dmcmcpar = as.double(dmcmc.par), imcmcpar = as.integer(imcmc.par), 
                               parmgsamp = double(nsamp * nparmg), acptsamp = double(nsamp * nblocks), lpsamp = double(nsamp), lprsamp = double(nsamp), 
                               lkappasamp = double(nsamp), phisamp = integer(nsamp), lpsisamp = double(nsamp), acpthpsamp = double(nsamp), 
                               acptlpsisamp = double(nsamp), lmglhsamp = double(nsamp), usamp = double(nsamp*n), fbase.choice = as.integer(fbase.choice)))
  if(verbose) cat("elapsed time:", round(tm.c[3]), "seconds\n")
  
  oo$x <- x; oo$y <- y; oo$phi.grid <- phi.grid; oo$location <- location
  oo$distance <- distance; oo$kappa <- kappa; oo$kernel <- kernel.choice
  oo$gridmats <- gridmats; oo$prox <- prox.grid; oo$reg.ix <- reg.ix;
  oo$runtime <- tm.c[3]
  class(oo) <- "qrjointtp"
  return(oo)
}
