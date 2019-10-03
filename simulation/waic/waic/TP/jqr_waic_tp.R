library(qrjoint)
waic.jqr <- matrix(NA, nrow = 100, ncol = 2)
for(i in 1:100){
  filename <- paste('jqrfit', i+600, '.RData', sep = '')
  load(file = filename)
  
  waic.jqr[i,] <- summary(jqr.fit, burn.perc = 0.5, ntrace = 5e2, plot.dev = FALSE, more.details = FALSE)$waic
  print(i)
}