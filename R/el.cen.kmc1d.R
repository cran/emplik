el.cen.kmc1d <- function(x,d,fun,mu, tol = .Machine$double.eps^0.5, step=0.001, ...) {
   xvec <- as.vector(x)
   nn <- length(xvec)
   if(nn <= 1) stop ("Need more observations")
   if(length(d)!=nn) stop("length of x and d must agree")
   if(any((d!=0)&(d!=1)))
     stop("d must be 0(right-censored) or 1(uncensored)")
   if(!is.numeric(xvec)) stop("x must be numeric")
   if(length(mu)!=1) stop("the dim of constraint, mu, must be 1")
   
   xd0 <- xvec[d==0]
   m <- length(xd0)

###############################################
#### do the computation in 2 different cases. #
###############################################
 
  if( (m <= 0) ) {          #### means no censor at all. 
        xd1 <- xvec
        funxd1 <- fun(xvec, ...)
        funNPMLE <- sum(funxd1)/nn    ###changed 10/2018
        logel00 <- - sum(log(nn))     ###changed 10/2018
        temp6 <- el.test.wt(funxd1, wt = rep(1,nn), mu=mu)   ### mu always 0? fixed.
        pnew <- temp6$prob
        lam <- temp6$lam
        logel <- sum(log(pnew))
  }
  
   if ((m > 0)) {
        tmp <- kmc.clean(Xtime = xvec, delta = d)
         sox <- tmp$kmc.time
         sod <- tmp$delta
         if (all(sod > 0) ) stop("the only censor is at front? (go to uncensor case after delet censor)")
         xd1 <- sox[sod == 1]
         xd0 <- sox[sod == 0]
         funxd1 <- fun(xd1, ...)
         gt <- fun(sox, ...) - mu        ##
        temp3 <- WKM(x = sox, d = sod, zc=1:length(sod), w = rep(1, length(sod)))
        logel00 <- temp3$logel
        funNPMLE <- sum(funxd1 * temp3$jump[temp3$jump > 0])
      
        m <- length(xd0)     #### may have changed
        n <- length(xd1)     ####length(pnew)
        k <- rep(NA, m) 
        for (i in 1:m) {
            k[i] <- 1 + n - sum(xd1 > xd0[i])
        }

       ff <- function(lambda, delta, gt) {
              temp4 <- omega.lambda(lambda=lambda, delta=delta, gt=gt)
              return(as.numeric(temp4$mea))
              }
        value0 <- ff(lambda=0.0, delta=sod, gt=gt)
        if( abs(value0) < tol) {lam <- 0}
            else{
                 if(value0 > 0) { stepLow <- -step
                                  temp2 <- uniroot(ff, lower=stepLow, upper=0, extendInt="upX", tol=tol, delta=sod, gt=gt)
                                 }
                 if(value0 < 0) { stepUp <- step
                                  temp2 <- uniroot(ff, lower=0, upper=stepUp, extendInt="upX", tol=tol, delta=sod, gt=gt)                 
                                 }
        ### suggested 10/17/2021 Mai Zhou
                 lam <- temp2$root
               }
       temp <- omega.lambda(lambda=lam, delta=sod, gt=gt)
       jumps <- temp$omega
       pnew <- jumps[jumps > 0]
       sur <- cumsumsurv(pnew)
       logel <- sum(log(pnew)) + sum(log(sur[k]))
    }
  
  #### need kmc.clean, WKM, omega.lambda, cumsumsurv

    tval <- 2 * (logel00 - logel)
    list(loglik = logel, times = xd1, prob = pnew, funMLE = funNPMLE, 
        lam = lam, `-2LLR` = tval, Pval = 1 - pchisq(tval, df = 1))
}
