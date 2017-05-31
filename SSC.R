## Copyright 2009, Joachim Werner, All Rights Reserved
## Address comments about this software to profwerner@t-online.de
##
## This software is made available AS IS. Warranty about the software, 
## its performance, or its conformity to any specification is NOT given 
## and NOT implied. The use of this software is permitted for academic 
## purposes only, provided its use is acknowledged.  For commercial use 
## the permission of the author is indispensable . 
##
## File       :  Signal summation conversion
## Calculates : "H"
## Updated    : Joachim Werner, Okt. 2009
##
## References:
## 1. Cannon, M.J., Percival, D.B., Caccia, D.C., Raymond, G.M. and Bassingthwaighte, J.B.,
##    Evaluating scaled windowed variance methods for estimating the Hurst coefficient of time series,
##    Physica A, 1997, 241, 606-626.
## 2. Eke, A., Herman, P., Bassingthwaighte, J.B., Raymond, G.M., Percival, D.B., Cannon, M., 
##    Balla, I. and Ikrenyi, C., 2000, Physiological time series: Distinguishing fractal 
##    noises from motions, Pfl�gers Archiv -- European Journal of Physiology, 439, 403-415.

SSC <- function(y, scale.min=4, scale.add=4, ...)   {
  #scale.min <- NULL            # minimale Intervalbreite
  #scale.add <- NULL            #### additive #### Vergr��erung der Intervalle (boxes)
  if (scale.min < 2) {stop("scale.min == 1 sinnfrei, deshalb STOP")}
  if (scale.add < scale.min) scale.add = scale.min 
  y <- cumsum(y - mean(y))      # Kumulierung der Zeitreihe
  n <- length(y)
  smax <- length(y)/2
  box <- scale.min              # box = variable Intervallbreite
  x <- 1:n                      # Generierung der x-Skala
  
  sd     <- rep(0, (n/scale.min))  # Fn
  scales <- rep(0, (n/scale.min))  # scales
  rep <- 0
  # Wiederholungen �ber variable Intervallbreiten
  repeat   {                         # Beginn unterschiedliche Intervallbreiten
    anz <- trunc(n/box)
    var <- 0
    for (i in 0:(anz-1))   {
      yr <- y[(1+i*box):(box+i*box)]
      
      bridge <- seq(yr[1], yr[box], length=box)    # a la R
      yr     <- yr - bridge                        # a la R
      yrm    <- yr - mean(yr)                                    
      var    <- var + sqrt(sum(yrm**2) / (box-1))
    }
    sdm <- var/anz                                 
    #       sdm <- sdm * n / (anz*box)  # Gewichtung, wenn n nicht Vielfaches von box!!! 
    #   cat(" rep", rep) 
    #   cat(" var", var)
    #   cat(" anz", anz)
    
    rep = rep + 1  
    
    sd[rep]     <- sdm
    scales[rep] <- box
    #   cat (" box", box)
    #   cat ("  sd[rep]", sd[rep], "\n")
    box <- box + scale.add
    if ( box >= smax) break
    #   if ( anz == 2) break
  }                        # Ende variable Intervallbreiten
  
  SD <-  log(sd[1:rep])
  SC <-  log(scales[1:rep])
  lmb <- lm(SD ~ SC)
  b <- coef(lmb)
  #cat (" H=", b[2])
  #plot(SC, SD, type="l")
  return(b[2])
}