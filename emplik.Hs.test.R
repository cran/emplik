
emplik.Hs.test <- function(y, d, gmat, theta, error=1e-8, maxit=20) {


#Some simulation results to confirm that for right censored data and k (k>1)
#hazard constraints, the limiting distribution is chi-square with k df.
#
#> set.seed(123)
#> for(i in 1:1000) result[i] <- simu(200)
#> for(i in 1001:5000) result[i] <- simu(200)
#> simu
#function(num) {
#
#xx <- rexp(num)
#cc <- 2*rexp(num)
#yy <- pmin(xx, cc)
#dd <- as.numeric(cc > xx)
#
#temp0 <- summary(survfit(Surv(yy,dd),se.fit=FALSE,conf.type="none" ) )

temp0 <- DnR(y,d) 
times1 <- temp0$time
k <- length(times1)
d1 <- temp0$n.event
R1 <- temp0$n.risk
if(d1[k]==R1[k]) {
                   d1<-d1[-k]
                   R1 <- R1[-k]
                   times1 <- times1[-k]
                   }

#col1 <- 0.5*times1*as.numeric(times1 < 1 )
#col2 <- exp(-times1)
#col3 <- as.numeric(times1 < 0.9)
#theta <- c(-0.25,-1,-0.9)
#
#gmat <- cbind(col1,col2,col3)
#
#iqpHazards(times1,d1,R1,gmat,theta)$ChiSq
#
#}
#
#
#temp <- DnR()
#
#iqpHazards <- function(times1, d1, R1, gmat, theta, error=1e-9, maxit=20) {

### This function computes the emp. lik. ratio for one sample right 
### censored data---with many hazard constraints. The input should be 
### in the style of summary(survfit())$times; n.events; n.risk. 
### which can be obtained by DnR() and Wdataclean2() too.
### The gmat is the matrix of g() values on the times1. 
### theta is a vector that the constraints are equal to:
### sum g(times1) log( 1- dH(times))=int g(t) dlog(1-F(t)) = theta
### Last changed  7/2001 by  Mai Zhou mai@ms.uky.edu 
### Also, d1/R1 must be < 1. If d1(last)=R1(last) then delete this obs.
### When the constraints are true,  -2EmpLikRatio should behave as Chisq (k).

w01 <- d1/R1
u01 <- log(1-w01)
mat0a <- (exp(-u01/2)-exp(u01/2))/sqrt(d1)
mat0b <- R1 - d1 - (d1*exp(u01))/(1-exp(u01))
mat0c <- as.matrix(gmat)
dd <- u01 %*% mat0c
k <- length(dd)
mat0d <- theta - dd

# value0 <- solve.QP(diag(mat0a),mat0b,mat0c,mat0d,meq=k,factorized=TRUE)
value0 <- solve3.QP(mat0a,mat0b,mat0c,mat0d,meq=k,factorized=TRUE)
wp <- u01 + value0$solution

###### iteration

m<-1
matd <- rep(0, k)
diff <- 10

while(diff > error && m < maxit) {
    mata <- (exp(-wp/2)-exp(wp/2))/sqrt(d1)
    matb <- R1-d1-(d1*exp(wp))/(1-exp(wp))
##  value <- solve.QP(diag(mata), matb, mat0c, matd, meq=k, factorized=TRUE)
 value <- solve3.QP(mata, matb, mat0c, matd, meq=k, factorized=TRUE)
    wp <- wp+value$solution
    diff <- sum(abs(value$solution))
    m<-m+1
    }
loglik2 <- sum( d1*log(1-exp(wp))+(R1-d1)*wp)
loglik1 <- sum(d1*log(w01)+(R1-d1)*u01)
list(weight = 1-exp(wp), ChiSq=2*(loglik1-loglik2))
}

#The Q-Q plot is in pp2.ps It agrees with chi square 3 df quite well.
