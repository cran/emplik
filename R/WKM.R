WKM <- function(x,d,w=rep(1, length(d)) ) {


temp <- Wdataclean2(x,d,w)
dd <- temp$dd
ww <- temp$weight
dd[length(dd)] <- 1

######why not use DnR?

allrisk <- rev(cumsum(rev(ww)))
survP <- cumprod( 1 -  (dd*ww)/allrisk )
jumps <- -diff( c(1, survP) )

logel <- sum(ww[dd==1]*log(jumps[dd==1])) + sum(ww[dd==0]*log(survP[dd==0]))

list(times=temp$value, jump=jumps, surv=survP, logel=logel)

}
################################################################
##This function compute a weighted Kaplan-Meier estimator
## x = times, d = censoring status, w = weights 
################################################################
