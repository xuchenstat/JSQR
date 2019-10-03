alqr <- function(tau, x, y, nphi = 20, alpha, update.alpha, sdlalpha = 10, decaypar, distance, kappa, kernel = c("se", "matern"), nsamp = 1e3, thin = 10, verbose = TRUE){
  
  n <- length(y)
  p <- ncol(x)
  niter <- nsamp * thin
  
  dimpars <- c(n, p, nphi, niter, thin, nsamp)

  cte <- 2^(1-kappa)/gamma(kappa)
  phi.lb <- uniroot(function(u) sapply(u, kernel, x = max(distance)/4, kappa = kappa, cte = cte) - 0.05, interval = c(1e-4, 10*max(distance)), tol = 1e-10)$root
  phi.ub <- 3*phi.lb
  phi.grid <- seq(phi.lb, phi.ub, length.out = nphi)
  cov.grid <- lapply(phi.grid, covmat, n = n, distance = distance, kernel = kernel, kappa = kappa, cte = cte)
  eigen.grid <- lapply(cov.grid, eigen)
  S <- sapply(eigen.grid, function(u) cbind(u$val, u$vec))
  
  phi <- floor(nphi/2)
  eta <- 1
  xi <- rep(1, n)
  
  sig.sq <- 1e4
  a0 <- 1e-3
  b0 <- 1e-3*sqrt(2/tau/(1-tau)+(1-2*tau)^2/tau^2/(1-tau)^2)
  
  init <- c(alpha, eta, xi)
  dhyper <- c(sig.sq, a0, b0, sdlalpha)

  tm.c <- system.time(oo <- .C("alp", tau = as.double(tau), x = as.double(x), y = as.double(y), S = as.double(S), dim = as.integer(dimpars),
                               phi = as.integer(phi), init = as.double(init), dhyper = as.double(dhyper), updateAlpha = as.integer(update.alpha),
                               betasamp = double(nsamp*p), etasamp = double(nsamp), xisamp = double(n*nsamp), phisamp = integer(nsamp), alphasamp = double(nsamp),
                               acptsamp = double(n+1)))
  oo$x <- x; oo$y <- y; oo$runtime <- tm.c[3]
  if(verbose) print(oo$runtime)
  return(oo)
}