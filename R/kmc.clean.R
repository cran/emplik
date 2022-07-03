kmc.clean <- function(Xtime, delta){
  ##TASK: 1. sort Time
  ##      2. delete the smallest censored obs. so that the first obs. is uncen! the last also uncen.

  n <- length(Xtime)
  if( all(delta == 0) ) stop('All obs. are censored?')
  dataOrder <- order(Xtime, -delta)
  kmc.time <- Xtime[dataOrder]
  delta <- delta[dataOrder]             ### changed 10/2018
  
  
  ####   tmp <- sort(kmc.time,index.return=TRUE)
  ####   kmc.time <- kmc.time[tmp$ix]
  ####   delta <- delta[tmp$ix]
  
  delta[n] <- 1
  FirstUnCenLocation <- which(delta==1)[1]
  if (FirstUnCenLocation==n) stop('Only one uncensored point.')
  if (FirstUnCenLocation!=1){
    delta <- delta[FirstUnCenLocation:n]
    kmc.time <- kmc.time[FirstUnCenLocation:n]
  }

  list(kmc.time=kmc.time, delta=delta)
}