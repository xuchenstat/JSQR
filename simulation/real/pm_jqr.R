pm <- read.csv('pm.csv') 
n <- nrow(pm)
n.train <- 270
n.test <- n - n.train

library(qrjoint)

for(datanum in 1:10){
  
  set.seed(datanum)
  
  train.ix <- sample(n, n.train)
  test.ix <- (1:n)[-train.ix]
  
  jqr.fit <- qrjoint(formula = pm ~ LCYPOP00 + LDISTA1 + NLCD_URB + LS25_E10 + ELEV_NED, data = pm[train.ix,],
                       nsamp = 2e3, thin = 10, fbase = 'logistic')
  
  save(jqr.fit, file = paste('jqrpm', datanum, '.RData', sep = ''))
}