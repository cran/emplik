omega.lambda <- function(lambda, delta, gt){
   n <- length(delta)    ###  delta is sorted according to kmc.time
   delta[n] <- 1         ###  Assume kmc.time already been kmc.clean-ed, so this is not needed.
   if( length(lambda) > 1 ) stop("lambda should be a scalar")
   lambda <- as.vector(lambda)   ### added 4/2020. not used in 1d, but for future.
   u.omega <- rep(0, n)   ##### all w = 0 to begin
   S <- rep(1, n)
  
   S.cen <- 0
   u.omega[1] <- 1/(n - lambda*gt[1])   
  
     for (k in 2:n){
       if(delta[k]==0){ S[k] <- S[k-1] - u.omega[k-1]
                        S.cen <- S.cen + 1/S[k] }
       else{ u.omega[k] <- 1/(n - lambda*gt[k] - S.cen)
             S[k] <- S[k-1] - u.omega[k-1] }
     }
  
###### temp <- .C('wlam', DDD=as.integer(delta), gti=as.numeric(gt), SSS=as.numeric(S), lam=as.numeric(lambda), LLL=as.integer(n))
######  
######    u.omega <- temp$gti
######    S <- temp$SSS
######  There are pure R version and C version.  
  
 list(Surv=S, omega=u.omega, mea= sum(gt*u.omega))
}




######################################
#### cumsumsurv <- function(x){
#### 	if(any(is.na(x))) stop('NaNs');    ## if (sum(is.na(x))>0) stop('NaNs');  3/2015 MZ
#### 	s=x;
#### 	.C('cumsumsurv',x=as.numeric(x),s=as.numeric(s),LLL=length(x))$s
#### 	}
###################################


