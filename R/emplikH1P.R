emplikH1P <- function(lambda, x, d, fung, CIforTheta=FALSE)
{
n <- length(x)
if( n <= 2 ) stop("Need more observations")
if( length(d) != n ) stop("length of x and d must agree")
if(any((d!=0)&(d!=1))) stop("d must be 0/1 for censor/not-censor")
if(!is.numeric(x)) stop("x must be numeric values -- observed times")
if(!is.numeric(lambda)) stop("lambda must be a numeric value -- tilting")

newdata <- Wdataclean2(x,d)
temp <- DnR(newdata$value, newdata$dd, newdata$weight)
time <- temp$times # only uncensored time? Yes.
risk <- temp$n.risk
jump <- (temp$n.event)/risk
funtime <- fung(time)
lam <- lambda

Nfuntime <- n * funtime
funh <- Nfuntime/risk # that is Zi
funtimeTjump <- funtime * jump
if(jump[length(jump)] >= 1) funh[length(jump)] <- 0 #for inthaz and weights

onepluslamh <- 1 + lam * funh ### this is 1 + lam Zi in Ref.
MeanHaz <- NA
if(CIforTheta) {
MeanHaz <- sum(funtimeTjump/onepluslamh)
}
## this should be used after the lambda interval is obtained,
## to turn it into theta interval.

weights <- jump/onepluslamh #need to change last jump to 1? NO. see above

loglik <- 2*(sum((temp$n.event)*log(onepluslamh)) - sum((temp$n.event)*(onepluslamh-1)/onepluslamh) )

#### sum(jump*(lam*sqrtfuntime)/onepluslamh)) ### allow ties 3/2020 ##  ##
##?is that right? YES see (3.2) in Ref. above. This is ALR, or Poisson LR.
list( "-2LLR"=loglik, ### logemlikv2=loglik2,
       lambda=lam, times=time, wts=weights, MeanHaz=MeanHaz )
}
