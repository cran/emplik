el.trun.test <- function(y, x, wt=rep(1,length(x)),
                  fun=function(x){x}, mu, maxit=20,error=1e-9) {
indi <- function(u,v){ as.numeric(u > v) }
uij <- outer(x,y,FUN="indi")
n <- length(x)
w0 <- wt/n
xmu <- fun(x) - mu
for(i in 1:maxit) {
     pvec0 <- as.vector( w0 %*% uij )
     nvec <- as.vector( rowsum( t(uij)*wt/pvec0, group=rep(1, n) ) )
     w0 <- wt/nvec
     w0 <- w0/sum(w0)
}
w <- w0
lamfun <- function(lam, nvec, xmu){sum(xmu/(nvec+lam*xmu))}
for(i in 1:maxit) {
       pvec <- as.vector( w %*% uij )
       nvec <- as.vector( rowsum( t(uij)*wt/pvec, group=rep(1, n) ) )
       BU <- 0.1*min(nvec)/max(abs(xmu))
       if(lamfun(0, nvec=nvec, xmu=xmu) == 0) lam <- 0
       else {
             if(lamfun(0, nvec=nvec, xmu=xmu) > 0) { lo <- 0
             up <- BU
             while(lamfun(up, nvec=nvec, xmu=xmu) > 0)
             up <- up + BU
             }
           else {up <- 0
                 lo <- -BU
                 while(lamfun(lo, nvec=nvec, xmu=xmu) < 0 ) lo <- lo - BU
                }
 lam <- uniroot(lamfun,lower=lo,upper=up, nvec=nvec,xmu=xmu,tol=1e-9)$root
            }
       w <- wt/(nvec + lam * xmu)
       w <- w/sum(w)
}
pvec <- as.vector( w %*% uij )
pvec0 <- as.vector( w0 %*% uij )
ELR <- sum(log(w0)-log(pvec0))-sum(log(w)-log(pvec))
return(list(NPMLE=w0, NPMLEmu=w, "-2LLR"=2*ELR) )
}
