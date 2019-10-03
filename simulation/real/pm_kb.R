pm <- read.csv('pm.csv') 
n <- nrow(pm)
n.train <- 270
n.test <- n - n.train


library(quantreg)

tau.grid <- c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99)


for(datanum in 1:10){
  
  set.seed(datanum)
  
  train.ix <- sample(n, n.train)
  test.ix <- (1:n)[-train.ix]
  
  kb.fit <- rq(formula = pm ~ LCYPOP00 + LDISTA1 + NLCD_URB + LS25_E10 + ELEV_NED, tau = tau.grid, data = pm[train.ix,])
  
  save(kb.fit, file = paste('kbpm', datanum, '.RData', sep = ''))
}
