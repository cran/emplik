WKM <- function(x,d,w=rep(1, length(d)) ) {

temp <- Wdataclean2(x,d,w)
dd <- temp$dd
ww <- temp$weight
dd[length(dd)] <- 1

allrisk <- rev(cumsum(rev(ww)))
survP <- cumprod( 1 -  (dd*ww)/allrisk )
jumps <- -diff( c(1, survP) )

list(times=temp$value, jump=jumps, surv=survP )

}

