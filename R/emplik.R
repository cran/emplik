####################################
####### el.test(), from Owen #######
####################################

el.test <- function( x, mu, lam, maxit=25, gradtol=1e-7, 
                        svdtol = 1e-9, itertrace=FALSE ){
x <- as.matrix(x)
n <- nrow(x)
p <- ncol(x)
mu <- as.vector(mu)
if( length(mu) !=p )
  stop("mu must have same dimension as observation vectors.")
if( n <= p )
  stop("Need more observations than length(mu) in el.test().")

z <- t( t(x) -mu )

#
#    Scale the problem, by a measure of the size of a 
# typical observation.  Add a tiny quantity to protect
# against dividing by zero in scaling.  Since z*lam is
# dimensionless, lam must be scaled inversely to z.
#
TINY <- sqrt( .Machine$double.xmin )
scale <- mean( abs(z) ) + TINY
z <- z/scale
if( !missing(lam) ){
  lam <- as.vector(lam)
  lam <- lam*scale
  if( logelr(z,rep(0,p),lam)>0 )lam <- rep(0,p)
}
if(  missing(lam)  )
  lam <- rep(0,p)
#
#     Take some precaution against users specifying
# tolerances too small.
#

if(  svdtol < TINY )svdtol <- TINY
if(  gradtol < TINY)gradtol <- TINY

#
#    Preset the weights for combining Newton and gradient
# steps at each of 16 inner iterations, starting with
# the Newton step and progressing towards shorter vectors
# in the gradient direction.  Most commonly only the Newton
# step is actually taken, though occasional step reductions
# do occur.
#

nwts <- c( 3^-c(0:3), rep(0,12) )
gwts <- 2^( -c(0:(length(nwts)-1)))
gwts <- (gwts^2 - nwts^2)^.5
gwts[12:16] <- gwts[12:16] * 10^-c(1:5)

#
#    Iterate, finding the Newton and gradient steps, and
# choosing a step that reduces the objective if possible.
#

nits <- 0
gsize <- gradtol + 1
while(  nits<maxit && gsize > gradtol  ){
  arg  <- 1 + z %*% lam
  wts1 <- as.vector( llogp(arg, 1/n) )
  wts2 <- as.vector( -llogpp(arg, 1/n) )^.5
  grad <- as.matrix( -z*wts1 )
  grad <- as.vector( apply( grad, 2, sum ) )
  gsize <- mean( abs(grad) )
  hess <- z*wts2
#                                   -1
#    The Newton step is -(hess'hess)    grad,
#  where the matrix hess is a sqrt of the Hessian.
#  Use svd on hess to get a stable solution.
#

  svdh <- svd( hess )
  if( min(svdh$d) < max(svdh$d)*svdtol )
    svdh$d <- svdh$d + max(svdh$d)*svdtol
  nstep <- svdh$v %*% (t(svdh$u)/svdh$d)
  nstep <- as.vector( nstep %*% matrix(wts1/wts2,n,1) )

  gstep <- -grad
  if(  sum(nstep^2) < sum(gstep^2) )
    gstep <- gstep*sum(nstep^2)^.5/sum(gstep^2)^.5
  ologelr <- -sum( llog(arg,1/n) )
  ninner <- 0
  for(  i in 1:length(nwts) ){
    nlogelr <- logelr( z,rep(0,p),lam+nwts[i]*nstep+gwts[i]*gstep )
    if( nlogelr < ologelr ){
      lam <- lam+nwts[i]*nstep+gwts[i]*gstep
      ninner <- i
      break
    }
  }
  nits <- nits+1
  if(  ninner==0  )nits <- maxit
  if( itertrace )
    print( c(lam, nlogelr, gsize, ninner) )
}

list( "-2LLR" = -2*nlogelr, Pval = 1-pchisq(-2*nlogelr, df=p),
     lambda = lam*scale, grad=grad*scale,
 hess=t(hess)%*%hess*scale^2, wts=wts1, nits=nits )
}


logelr <- function( x, mu, lam ){
x <- as.matrix(x)
n <- nrow(x)
p <- ncol(x)
if(  n <= p  )
  stop("Need more observations than variables in logelr.")
mu <- as.vector(mu)
if(  length(mu) != p  )
  stop("Length of mean doesn't match number of variables in logelr.")

z <- t( t(x) -mu )
arg <- 1 + z %*% lam
- sum( llog(arg,1/n) )
}

#
#    The function llog() is equal to the natural
#  logarithm on the interval from eps >0 to infinity.
#  Between -infinity and eps, llog() is a quadratic.
#  llogp() and llogpp() are the first two derivatives
#  of llog().  All three functions are continuous
#  across the "knot" at eps.
#
#    A variation with a second knot at a large value
#  M did not appear to work as well.
#
#    The cutoff point, eps, is usually 1/n, where n
#  is the number of observations.  Unless n is extraordinarily
#  large, dividing by eps is not expected to cause numerical
#  difficulty.
#

llog <- function( z, eps ){

ans <- z
lo <- (z<eps)
ans[ lo  ] <- log(eps) - 1.5 + 2*z[lo]/eps - 0.5*(z[lo]/eps)^2
ans[ !lo ] <- log( z[!lo] )
ans
}

llogp <- function( z, eps ){

ans <- z
lo <- (z<eps)
ans[ lo  ] <- 2.0/eps - z[lo]/eps^2
ans[ !lo ] <- 1/z[!lo]
ans
}

llogpp <- function( z, eps ){

ans <- z
lo <- (z<eps)
ans[ lo  ] <- -1.0/eps^2
ans[ !lo ] <- -1.0/z[!lo]^2
ans
}

###################################
#######  emplikH1.test() ##########
###################################

emplikH1.test <- function(x, d, theta, fun, 
	              tola = .Machine$double.eps^.25)
{
n <- length(x)
if( n <= 2 ) stop("Need more observations")
if( length(d) != n ) stop("length of x and d must agree")
if(any((d!=0)&(d!=1))) stop("d must be 0/1's for censor/not-censor")
if(!is.numeric(x)) stop("x must be numeric values --- observed times")

#temp<-summary(survfit(Surv(x,d),se.fit=F,type="fleming",conf.type="none"))
#
newdata <- Wdataclean2(x,d)
temp <- DnR(newdata$value, newdata$dd, newdata$weight)

time <- temp$time         # only uncensored time?  Yes. 
risk <- temp$n.risk
jump <- (temp$n.event)/risk

funtime <- fun(time)
funh <- (n/risk) * funtime    # that is Zi 
funtimeTjump <- funtime * jump 

if(jump[length(jump)] >= 1) funh[length(jump)] <- 0  #for inthaz and weights

inthaz <- function(x, ftj, fh, thet){ sum(ftj/(1 + x * fh)) - thet } 

diff <- inthaz(0, funtimeTjump, funh, theta)

if( diff == 0 ) { lam <- 0 } else {
    step <- 0.2/sqrt(n) 
    if(abs(diff) > 6*log(n)*step )
    stop("given theta value is too far away from theta0")

    mini<-0
    maxi<-0 
    if(diff > 0) { 
    maxi <- step 
    while(inthaz(maxi, funtimeTjump, funh, theta) > 0  && maxi < 50*log(n)*step) 
    maxi <- maxi+step 
    } 
    else { 
    mini <- -step 
    while(inthaz(mini, funtimeTjump, funh, theta) < 0 && mini > - 50*log(n)*step) 
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
#?is that right? YES  see (3.2) in Ref. above. This ALR, or Poisson LR.

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

list( logemlik=loglik,  ### logemlikv2=loglik2, 
      lambda=lam, times=time, wts=weights, 
      nits=temp2$nf, message=temp2$message )
}

Wdataclean2 <- function (z, d, wt = rep(1,length(z)) ) 
{  niceorder <- order(z,-d)
   sortedz <- z[niceorder]
   sortedd <- d[niceorder]
   sortedw <- wt[niceorder]

   n <- length(sortedd)
   y1 <- sortedz[-1] != sortedz[-n]
   y2 <- sortedd[-1] != sortedd[-n]
   y <- y1 | y2

   ind <- c(which(y | is.na(y)), n)

   csumw <- cumsum(sortedw)

   list( value = sortedz[ind], dd = sortedd[ind],
         weight = diff(c(0, csumw[ind])) )
}

DnR <- function(x, d, w)
{
# inputs should be from  Wdataclean2()

allrisk <- rev(cumsum(rev(w)))
posi <- d == 1
uncenx <- x[posi]
uncenw <- w[posi]
uncenR <- allrisk[posi]

list( time = uncenx, n.risk = uncenR, n.event = uncenw )
}


###############################################
############# emplikH2() ######################
###############################################

emplikH2.test <- function(x, d, K, fun, 
	                tola = .Machine$double.eps^.25,...)
{
if(!is.numeric(x)) stop("x must be numeric values -- observed times")
n <- length(x) 
if( n <= 2 ) stop("Need more observations than two")
if( length(d) != n ) stop("length of x and d must agree") 
if(any((d!=0)&(d!=1))) stop("d must be 0/1's for censor/not-censor")

#temp <- summary(survfit(Surv(x,d),se.fit=F,type="fleming",conf.type="none"))
#
newdata <- Wdataclean2(x,d)
temp <- DnR(newdata$value, newdata$dd, newdata$weight)


Dtime <- temp$time         # only uncensored time?  Yes. 
risk <- temp$n.risk 
jump <- (temp$n.event)/risk

funtime <- fun(Dtime,...)
funh <- (n/risk) * funtime                      # that is Zi  
funtimeTjump <- funtime * jump 

if(jump[length(jump)] >= 1) funh[length(jump)] <- 0  #for inthaz and weights

inthaz <- function(x, ftj, fh, Konst){ sum(ftj/(1 + x * fh)) - Konst } 

diff <- inthaz(0, funtimeTjump, funh, K)

if( diff == 0 ) { lam <- 0 } else {
    step <- 0.2/sqrt(n) 
    if(abs(diff) > 99*log(n)*step )        ##why 99*log(n)? no reason, you 
    stop("given theta value is too far away from theta0") # need something. 

    mini<-0
    maxi<-0 
    if(diff > 0) { 
    maxi <- step 
    while(inthaz(maxi, funtimeTjump, funh, K) > 0  && maxi < 50*log(n)*step) 
    maxi <- maxi+step 
    } 
    else { 
    mini <- -step 
    while(inthaz(mini, funtimeTjump, funh, K) < 0 && mini > - 50*log(n)*step) 
    mini <- mini - step 
    } 

    if(inthaz(mini,funtimeTjump,funh,K)*inthaz(maxi,funtimeTjump,funh,K)>0)
    stop("given theta is too far away from theta0")

    temp2 <- uniroot(inthaz,c(mini,maxi), tol = tola, 
                  ftj=funtimeTjump, fh=funh, Konst=K)  
    lam <- temp2$root 
} 

onepluslamh<- 1 + lam * funh   # this is 1 + lam Zi in Ref. 

weights <- jump/onepluslamh  #need to change last jump to 1? NO. see above

loglik <- 2*(sum(log(onepluslamh)) - sum((onepluslamh-1)/onepluslamh) )
#?is that right? YES see (3.2) in Ref. above. This is ALR, or Poisson LR.

#last <- length(jump)      #to compute loglik2, we need to drop last jump
#if (jump[last] == 1) {
#                     risk1 <- risk[-last]
#                     jump1 <- jump[-last]
#                     weights1 <- weights[-last]
#                     } else {
#                            risk1 <- risk
#                            jump1 <- jump
#                            weights1 <- weights
#                            }
#
#loglik2 <- 2*( sum(log(onepluslamh)) + 
#          sum( (risk1 -1)*log((1-jump1)/(1- weights1) ) )  ) 
# this version of LR seems to give negative value sometime???

list( "-2logemLR"=loglik,  ### drop this output "-2logemLRv2"=loglik2, 
      lambda=lam, times=Dtime, wts=weights, 
      nits=temp2$nf, message=temp2$message )
}
# what should be the fun() and K if I want to perform a (1-sample) 
# log-rank test?
# fun3 <- function(t1, z1) { psum( t( outer(z1, t1, FUN=">=") ) ) } 
# this is similar to the function in LogRank2() function. Need psum/2/3.
# And K = int R(t) dH(t)  = sum( H(z1) ) For example if H() is 
# exponential(0.3) then H(t) = 0.3*t, i.e. K <- sum(0.3* z1) 
# so finally a call may look like
#
# Assume z1 and d1 are the inputs:
# emlik2(z1, d1, sum(0.3* z1), 
#   fun3 <- function(t1,z){psum(t(outer(z,t1,FUN=">=") ) ) }, z=z1)
#
# Now use z1<-c(1,2,3,4,5) and d1<-c(1,1,0,1,1) we get
# emlik2(z1, d1, sum(0.25* z1),
#   fun3 <- function(t1,z){psum(t(outer(z,t1,FUN=">=") ) ) }, z=z1)
#
# with outputs that include this (and more)
# $ "-2logemLR":
# [1] 0.02204689
#This tests if the (censored) obs. z1 is from exp(0.25)


########################################################
############ emplikdisc.test() #########################
########################################################

emplikdisc.test <- function(x, d, K, fun, 
	                     tola=.Machine$double.eps^.25, theta)
{
n <- length(x) 
if(n <= 2) stop("Need more observations")
if(length(d) != n ) stop("length of x and d must agree") 
if(any((d!=0)&(d!=1))) stop("d must be 0/1's for censor/not-censor")
if(!is.numeric(x)) stop("x must be numeric values --- observed times")

#temp<-summary(survfit(Surv(x,d),se.fit=F,type="fleming",conf.type="none"))
#
newdata <- Wdataclean2(x,d)
temp <- DnR(newdata$value, newdata$dd, newdata$weight)

otime <- temp$time         # only uncensored time?  Yes. 
orisk <- temp$n.risk
odti <- temp$n.event

###if the last jump is of size 1, we need to drop last jump from computation
last <- length(orisk) 
if (orisk[last] == odti[last]) {
                     otime <- otime[-last] 
                     orisk <- orisk[-last]
                     odti  <- odti[-last]
                     }
######## compute the function g(ti, theta) 
gti <- fun(otime,theta) 

### the constrain function. To be solved in equation later.

constr <- function(x, Konst, gti, rti, dti, n) { 
                  rtiLgti <- rti + x*n*gti
                  OneminusdH <- (rtiLgti - dti)/rtiLgti
                  if( any(OneminusdH <= 0) ) stop(" too far away ")
                  sum(gti*log(OneminusdH)) -  Konst } 

##############################################################

differ <- constr(0, Konst=K, gti=gti, rti=orisk, dti=odti, n=n)

if( abs(differ) < tola ) { lam <- 0 } else {
    step <- 0.2/sqrt(n) 
    if(abs(differ) > 200*log(n)*step )   #Why 60 ? 
    stop("given theta value is too far away from theta0")

    mini<-0
    maxi<-0            
######### assume the constrain function is increasing in lam (=x) 
    if(differ > 0) { 
    mini <- -step 
    while(constr(mini, Konst=K, gti=gti, rti=orisk, dti=odti, n=n) > 0
 	          && mini > -200*log(n)*step ) 
    mini <- mini - step 
    } 
    else { 
    maxi <- step 
    while(constr(maxi, Konst=K, gti=gti, rti=orisk, dti=odti, n=n) < 0
                  &&  maxi < 200*log(n)*step ) 
    maxi <- maxi+step 
    }

    if(constr(mini, Konst=K, gti=gti, rti=orisk, dti=odti, n=n)*constr(maxi, 
                Konst=K, gti=gti, rti=orisk, dti=odti, n=n) > 0 )
    stop("given theta/K is/are too far away from theta0/K0")

# Now we solve the equation to get lambda, to satisfy the constraint of Ho

    temp2 <- uniroot(constr,c(mini,maxi), tol = tola, 
                  Konst=K, gti=gti, rti=orisk, dti=odti, n=n)  
    lam <- temp2$root 
}
####################################################################
rPlgti <- orisk + n*lam*gti

loglik <- 2*sum(odti*log(rPlgti/orisk) +
           (orisk-odti)*log(((orisk-odti)*rPlgti)/(orisk*(rPlgti-odti)) ) )

#?is that right? YES the -2log lik ratio. 
# Notice the output time and jumps has less the last point.
list("discrete.-2logemlikRatio"=loglik, lambda=lam, times=otime,
                jumps=odti/rPlgti)
}

################################################
####### solve3.QP, WKM, and el.cen.test ########
################################################

solve3.QP <- function(D, d, A, b, meq, factorized=FALSE) {
#### This code works for QP problem: min 1/2x'Dx-d'x
#### where the matrix D is diagonal and the constraints
#### are all equalities, i.e. t(A)x=b.
#### Inputs:
#### D should be a vector of length n, this means the matrix diag(D), but
#### if factorized=TRUE, D actually is diag(D)^(-1/2).
#### d is a vector of length n
#### A is a matrix of n x p
#### b is a vector with length p.  Finally meq = integer p.
#### The input meq are here for the compatibility with solve.QP in R package
#### quadprog. Written by Mai Zhou (mai@ms.uky.edu) Jan.30, 2001
D <- as.vector(D)
if(length(b)!=meq) stop("length of constraints not matched")
if(length(D)!=length(d)) stop("dimention of D and d not match")
if(dim(A)[1]!=length(D)) stop("dimention of D and A not match")
if(dim(A)[2]!=meq) stop("dimention of A not match with meq")

         if(!factorized) { D <- 1/sqrt(D) }
         QRout <- qr(D*A)
         temp <- rep(0,meq)
         if(any(b!=0)) {temp<-forwardsolve(t(qr.R(QRout)),b)}
         temp2 <- temp - t(qr.Q(QRout)) %*% (D * d)
         eta <- backsolve(qr.R(QRout), temp2)
         sol <- D^2 * (d + A %*% eta )
list( solution=sol )
}

el.cen.test <- function(x,d,fun=function(x){x},mu,error=1e-8,maxit=20)
{
   xvec <- as.vector(x)
   n <- length(xvec)
   if(n <= 2) stop ("Need more observation")
   if(length(d)!=n) stop("length of x and d must agree")
   if(any((d!=0)&(d!=1)))
      stop ("d must be 0(right-censored) or 1(uncensored)")
   if(!is.numeric(xvec)) stop("x must be numeric")
   if(length(mu)!=1) stop("check the dim of constraint mu")

   temp <- Wdataclean2(xvec,d)
   dd <- temp$dd
   dd[length(dd)] <- 1

   if(all(dd==1)) stop("there is no censoring, please use el.test()")

   xx <- temp$value
   n <- length(xx) 
   ww <- temp$weight
   w0 <- WKM(xx, dd, ww)$jump
   uncenw0 <- w0[dd==1]
   funxx <- fun(xx)

   if((mu>max(funxx))|(mu<min(funxx))) stop("check the value of mu/fun")

   xbar <- sum(funxx[dd==1] * uncenw0)

    #********* begin initial calculation******************
    # get vector dvec which is the first derivative vector

     dvec01 <- uncenw0
     rk <- 1:n            ######## rank(sortx) = 1:n  yes!
     cenrk <- rk[dd==0]
     mm <- length(cenrk)
     dvec02 <- rep(0,mm)
     for(j in 1:mm)  dvec02[j] <- sum(w0[cenrk[j]:n])
     dvec00 <- rep(0,n)
     dvec00[dd==1] <- dvec01
     dvec00[dd==0] <- dvec02
     dvec0 <- ww/dvec00

     # get matix Dmat which is Decompition of 2nd derivative matrix.
     # Dmat0 <- diag(1/dvec0)
     Dmat0 <- dvec00/sqrt(ww)

     # get constraint matrix Amat
     mat <- matrix(rep(dd,mm),ncol=mm, nrow=n)
     for(i in 1:mm)
     {
        mat[1:cenrk[i],i] <- 0
        mat[cenrk[i],i] <- -1
     }

    Amat <- as.matrix(cbind(dd, funxx * dd, mat))

     # get constraint vector bvec
     bvec0 <- c(0,as.vector(mu-xbar),rep(0,mm))

     # Use solve3.QP to maximize the loglikelihood function
     value0<-solve3.QP(Dmat0,dvec0,Amat,bvec0,meq=mm+2,factorized=TRUE)

     w <- dvec00 + value0$solution

     if(any(w<=0))
     stop("There is no probability satisfying the constraints")

     #**********end initial calculation **********************
     #**********begin iteration ******************************
     # update vector Dmat, dvec and bvec after initial calculation
     bvec <- rep(0,mm+2)
     diff <- 10
     m <- 0
     while( (diff>error) & (m<maxit) )
      {
         dvec <- ww/w
         # get matix Dmat
         #  Dmat <- diag(w)
         Dmat <- w/sqrt(ww)
         value0 <- solve3.QP(Dmat,dvec,Amat,bvec,meq=mm+2,factorized=TRUE)
         w <- w + value0$solution
         #diff <- sum(value0$solution^2)
         diff <- sum(abs( value0$solution) )
         m <- m+1
       }
     #**********end iteration ******************************

     lik00 <- sum(ww*log(dvec00))

   list("-2LLR"=2*(lik00 - sum(ww*log(w))), weights=w[dd==1],
                    xtime=xx[dd==1], iteration=m, error=diff)
}

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

