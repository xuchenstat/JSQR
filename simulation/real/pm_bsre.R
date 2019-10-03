pm <- read.csv('pm.csv') 
coords <- pm[,c('x','y')]

library(spBayes)

distance <- as.vector(dist(coords, method = 'euclidean'))

matern <- function(x, phi, kappa, cte){
  uphi <- x/phi;
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


kappa = 2
cte <- 2^(1-kappa)/gamma(kappa)
phi.lb <- uniroot(function(u) sapply(u, matern, x = max(distance)/4, kappa = kappa, cte = cte) - 0.05, interval = c(1e-4, 10*max(distance)), tol = 1e-10)$root
phi.ub <- 3*phi.lb

priors <- list('beta.Flat', 'tau.sq.IG'=c(1,1), 'sigma.sq.IG' = c(1,1), 'nu.Unif' = c(2-1e-10,2+1e-10), 
               'phi.Unif' = 1/c(phi.ub, phi.lb))
starting <- list('tau.sq' = 1, 'sigma.sq' = 1, 'phi' = mean(1/c(phi.lb,phi.ub)), 'nu' = 2)

tuning <- list('tau.sq' = 0.01, 'sigma.sq' = 0.01, 'phi'= 1.5, 'nu'=0)

set.seed(nrow(pm))
bsre.fit <- spLM(formula = pm ~ LCYPOP00 + LDISTA1 + NLCD_URB + LS25_E10 + ELEV_NED, data = pm, coords = as.matrix(coords), 
                 cov.model = 'matern', starting = starting, priors = priors, tuning = tuning, n.samples = 2e4, n.report = 2e3)


save(bsre.fit, file = 'bsrepm.RData')