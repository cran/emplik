emplikH1B <- function(lambda, x, d, fung, CIforTheta=FALSE)
{
n <- length(x)
if(n <= 2) stop("Need more observations")
if(length(d) != n ) stop("length of x and d must agree")
if(any((d!=0)&(d!=1))) stop("d must be 0/1 for censor/not-censor")
if(!is.numeric(x)) stop("x must be numeric values -- observed times")
if(!is.numeric(lambda)) stop("lambda must be a numeric value -- tilt")

newdata <- Wdataclean2(x,d)
newdata$dd[length(newdata$dd)] <- 1    #### added July 2022 MZ

temp <- DnR(newdata$value, newdata$dd, newdata$weight)
otime <- temp$times # only uncensored time? Yes.
orisk <- temp$n.risk
odti <- temp$n.event
###if the last jump is of size 1, we need to drop last jump from computation
last <- length(orisk)
if (orisk[last] == odti[last]) {
otime <- otime[-last]
orisk <- orisk[-last]
odti <- odti[-last]
}

######## compute the function g(ti)
gti <- fung(otime)
lam <- lambda
########### the constrain function.
constr <- function(lam, gti, rti, dti, n) {
                   rtiLgti <- rti + lam*n*gti
                   OneminusdH <- (rtiLgti - dti)/rtiLgti
                   if( any(OneminusdH <= 0) ) stop("estimator not well defined")
                   return(-sum(gti*log(OneminusdH))) 
}

IntHaz <- NA
if(CIforTheta) {
IntHaz <- constr(lam=lam, gti=gti, rti=orisk, dti=odti, n=n)
}
## Kcenter <- sum(gti * log(1- odti/orisk)) ## not used??
##############################################################
rPlgti <- orisk + n*lam*gti
loglik <- 2*sum(odti*log(rPlgti/orisk) + (orisk-odti)*log(((orisk-odti)*rPlgti)/(orisk*(rPlgti-odti)) ) )

## Notice the output time and jumps has less the last point.
list("-2LLR"=loglik, lambda=lam, times=otime,
jumps=odti/rPlgti, IntHaz=IntHaz)
}