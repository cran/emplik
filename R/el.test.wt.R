el.test.wt <- function(x, wt, mu) {
xmu <- x-mu
allw <- sum(wt)
BU <- 0.1*allw/max(abs(xmu))

lamfun <- function(lam,xmu,wt,allw) { sum(wt*xmu/(allw+lam*xmu)) }

if(lamfun(0,xmu,wt,allw) == 0) lam0 <- 0 
else {
 if( lamfun(0,xmu,wt,allw) > 0 ) {lo=0
                                up <- BU
                                while(lamfun(up,xmu,wt,allw)>0)
                                     up <- up + BU
                                 }
 else {up=0
      lo= - BU
      while(lamfun(lo,xmu,wt,allw) < 0 )
           lo <- lo - BU
     }
 lam0 <- uniroot(lamfun,lower=lo,upper=up,tol=1e-9,xmu=xmu,wt=wt,allw=allw)$root
}
pi <- wt/(allw + lam0*xmu)
list(x=x, wt=wt, prob=pi)
}
