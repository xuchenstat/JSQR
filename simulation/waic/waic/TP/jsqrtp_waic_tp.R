dyn.load('/jsqr/jsqrtp/jsqrtp.so')
source('/jsqr/jsqrtp/utility.R')
source('/jsqr/jsqrtp/jsqrtp.R')

waic.jsqrtp <- matrix(NA, nrow = 100, ncol = 2)

for(i in 1:100){
  filename <- paste('jsqrtpfit', i+600, '.RData', sep = '')
  load(file = filename)
  
  waic.jsqrtp[i,] <- waictp(jsqr.fit, burn.perc = 0.5, nmc = 500)
  print(i)
}