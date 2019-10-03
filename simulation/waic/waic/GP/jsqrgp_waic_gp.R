dyn.load('/jsqr/jsqrgp/jsqrgp.so')
source('/jsqr/jsqrgp/utility.R')
source('/jsqr/jsqrgp/jsqrgp.R')

waic.jsqrgp <- matrix(NA, nrow = 100, ncol = 2)

for(i in 1:100){
  filename <- paste('jsqrgpfit', i+500, '.RData', sep = '')
  load(file = filename)
  
  waic.jsqrgp[i,] <- waicgp(jsqr.fit, burn.perc = 0.5, nmc = 500)
  print(i)
}
