kmc.clean <- function(kmc.time, delta){
  ##TASK: 1 sort Time
  ##2 the first is uncen! the last also uncen.

  n <- length(kmc.time)
  dataOrder <- order(kmc.time, -delta)
  kmc.time <- kmc.time[dataOrder]
  delta <- delta[dataOrder]             ### changed 10/2018
  
  ####   tmp <- sort(kmc.time,index.return=TRUE)
  ####   kmc.time <- kmc.time[tmp$ix]
  ####   delta <- delta[tmp$ix]
  
  delta[n] <- 1
  FirstUnCenLocation <- which(delta==1)[1]
  if (FirstUnCenLocation==n) {stop('Only one uncensored point.')}
  if (FirstUnCenLocation!=1){
    delta <- delta[FirstUnCenLocation:n]
    kmc.time <- kmc.time[FirstUnCenLocation:n]
  }

  return (list(kmc.time=kmc.time,delta=delta))
}