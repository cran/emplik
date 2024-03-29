###################################
#######  emplikH1.test() ##########
###################################

emplikH1.test <- function(x, d, y= -Inf, theta, fun, 
	              tola = .Machine$double.eps^.5)
{
n <- length(x)
if( n <= 2 ) stop("Need more observations")
if( length(d) != n ) stop("length of x and d must agree")
if(any((d!=0)&(d!=1))) stop("d must be 0/1's for censor/not-censor")
if(!is.numeric(x)) stop("x must be numeric values --- observed times")

#
#temp<-summary(survfit(Surv(x,d),se.fit=F,type="fleming",conf.type="none"))
#
newdata <- Wdataclean2(x,d)
temp <- DnR(newdata$value, newdata$dd, newdata$weight, y=y)

time <- temp$times         # only uncensored time?  Yes. 
risk <- temp$n.risk
jump <- (temp$n.event)/risk

funtime <- fun(time)
funh <- sqrt(n) * funtime/risk    # that is Zi/sqrt(n) 
funtimeTjump <- funtime * jump 

if(jump[length(jump)] >= 1) funh[length(jump)] <- 0  #for inthaz and weights

inthaz <- function(x, ftj, fh, thet){ sum(ftj/(1 + x * fh)) - thet } 

diff <- inthaz(0, funtimeTjump, funh, theta)

if( diff == 0 ) { lam <- 0 } else {
    step <- 0.2/sqrt(n) 
    ####  if(abs(diff) > 50*log(n)*step )
    ####  stop("given theta value is too far away from theta0")

    mini<-0
    maxi<-0
    if(diff > 0) {
    maxi <- step
    while(inthaz(maxi, funtimeTjump, funh, theta) > 0 && maxi < 300*log(n)*step) 
    maxi <- maxi+step 
    }                    ##### maximum steps increased from 50 to 300, two places. 9/2022 MZ
    else { 
    mini <- -step 
    while(inthaz(mini, funtimeTjump, funh, theta) < 0 && mini > - 300*log(n)*step) 
    mini <- mini - step 
    } 

    if(inthaz(mini, funtimeTjump, funh, theta)*inthaz(maxi, funtimeTjump, funh, theta) > 0 )
    stop("given theta is too far away from theta0")

    temp2 <- uniroot(inthaz,c(mini,maxi), tol = tola, 
                  ftj=funtimeTjump, fh=funh, thet=theta)  
    lam <- temp2$root
}

onepluslamh<- 1 + lam * funh   ### this is 1 + lam Zi in Ref. 

weights <- jump/onepluslamh  #need to change last jump to 1? NO. see above

loglik <- 2*(sum(log(onepluslamh)) - sum((onepluslamh-1)/onepluslamh) ) 
#?is that right? YES  see (3.2) in Ref. above. This is ALR, or Poisson LR.
# log() function may generate NaN.
#last <- length(jump)    ## to compute loglik2, we need to drop last jump
#if (jump[last] == 1) {
#                     risk1 <- risk[-last]
#                     jump1 <- jump[-last]
#                     weights1 <- weights[-last]
#                     } else { 
#                            risk1 <- risk
#                            jump1 <- jump
#                            weights1 <- weights
#                            }
#loglik2 <- 2*( sum(log(onepluslamh)) + 
#          sum( (risk1 -1)*log((1-jump1)/(1- weights1) ) )  ) 
##? this likelihood seems have negative values sometimes???

list( "-2LLR"=loglik,  ### logemlikv2=loglik2, 
      lambda=lam/sqrt(n), times=time, wts=weights, 
      nits=temp2$nf, message=temp2$message )
}

