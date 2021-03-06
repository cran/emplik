\name{el.cen.EM2}
\alias{el.cen.EM2}
\title{Empirical likelihood ratio test for a vector of means 
with right, left or doubly censored data, by EM algorithm}
\usage{
el.cen.EM2(x,d,xc=1:length(x),fun,mu,maxit=50,error=1e-9,...)
}
\description{
This function is similar to \code{el.cen.EM()}, but for multiple constraints.
In the input there is a vector of observations 
\eqn{x = (x_1, \cdots , x_n)} and a 
function \code{fun}. The function \code{fun} should return the 
(n by k) matrix
\deqn{
          ( f_1(x), f_2(x), \cdots, f_k (x) ) . 
}

Also, the ordering of the observations, when consider censoring or 
redistributing-to-the-right, 
is according to the value of \code{x}, not \code{fun(x)}. 
So the probability distribution is for values \code{x}.
This program uses EM algorithm to maximize 
(wrt \eqn{p_i}) empirical
log likelihood function for right, left or doubly censored data with 
the MEAN constraint:
\deqn{ j = 1,2, \cdots ,k ~~~~ 
   \sum_{d_i=1} p_i f_j(x_i) = \int f_j(t) dF(t) = \mu_j ~. }
Where \eqn{p_i = \Delta F(x_i)} is a probability,
\eqn{d_i} is the censoring indicator, 1(uncensored), 0(right censored),
2(left censored). 
It also returns those \eqn{p_i}. 
The log likelihood function is defined as
\deqn{ \sum_{d_i=1} \log \Delta F(x_i)  + \sum_{d_i=2} \log F(x_i) 
     + \sum_{d_i=0} \log [ 1-F(x_i)] ~.}
}
\arguments{
    \item{x}{a vector containing the observed survival times.}
    \item{d}{a vector containing the censoring indicators, 
           1-uncensored; 0-right censored; 2-left censored.}
    \item{xc}{an optional vector of collapsing control values. 
              If xc[i] xc[j] have different values then 
              (x[i], d[i]), (x[j], d[j]) will not merge into one 
              observation with weight two, even 
              if they are identical. Default is not to merge.}
    \item{fun}{a left continuous (weight) function that returns a matrix. 
         The columns (=k) of the matrix is used to calculate
         the means and will be tested in \eqn{H_0}.
         \code{fun(t)} must be able to take a vector input \code{t}.} 
    \item{mu}{a vector of length k. Used in the constraint, 
                    as the mean of \eqn{f(X)}.}
    \item{maxit}{an optional integer, used to control maximum number of
             iterations. }
    \item{error}{an optional positive real number specifying the tolerance of
       iteration error. This is the bound of the
       \eqn{L_1} norm of the difference of two successive weights.}
    \item{...}{additional inputs to pass to \code{fun()}.}
}
\value{
    A list with the following components:
    \item{loglik}{the maximized empirical log likelihood under the constraints.}
    \item{times}{locations of CDF that have positive mass.}
    \item{prob}{the jump size of CDF at those locations.}
    \item{"-2LLR"}{If available, it is Minus two times the 
                   Empirical Log Likelihood Ratio.
                   Should be approx. chi-square distributed under Ho.}
    \item{Pval}{If available, the P-value of the test, 
               using chi-square approximation.}
    \item{lam}{the Lagrange multiplier in the final EM step. (the M-step)}
}
\details{

This implementation is all in R and have several for-loops in it. 
A faster version would use C to do the for-loop part.
(but this version is easier to port to Splus, and seems faster enough). 

We return the log likelihood all the time. Sometimes, (for right censored
and no censor case) we also return the -2 log likelihood ratio.
In other cases, you have to plot a curve with many values of the 
parameter, mu, to
find out where the log likelihood becomes maximum.
And from there you can get -2 log likelihood ratio between
the maximum location and your current parameter in Ho.

In order to get a proper distribution as NPMLE, we automatically
change the \eqn{d} for the largest observation to 1
(even if it is right censored), similar for the left censored, 
smallest observation.
\eqn{\mu} is a given constant vector. 
When the given constants \eqn{\mu} is too far
away from the NPMLE, there will be no distribution
satisfy the constraint.
In this case the computation will stop.
The -2 Log empirical likelihood ratio
should be infinite. 

The constant vector \code{mu} must be inside 
\eqn{( \min f(x_i) , \max f(x_i) ) }
for the computation to continue. 
It is always true that the NPMLE values are feasible. So when the
computation stops, try move the \code{mu} closer
to the NPMLE --- 
\deqn{ \hat \mu _j = \sum_{d_i=1} p_i^0 f_j(x_i) } 
where \eqn{p_i^0} taken to be the jumps of the NPMLE of CDF. 
Or use a different \code{fun}. 

Difference to the function \code{el.cen.EM}: due to the introduction of
input \code{xc} here in this function, the output \code{loglik} may be different
compared to the function \code{el.cen.EM}
due to not collapsing of duplicated input survival values.
The \code{-2LLR} should be the same from both functions.

}
\author{ Mai Zhou }
\references{
    Zhou, M. (2005). Empirical likelihood ratio with arbitrary censored/truncated data by EM algorithm. 
 \emph{Journal of Computational and Graphical Statistics}, 643-656.

    Zhou, M. (2002). 
        Computing censored empirical likelihood ratio 
        by EM algorithm. 
    \emph{Tech Report, Univ. of Kentucky, Dept of Statistics}

}
\examples{
## censored regression with one right censored observation.
## we check the estimation equation, with the MLE inside myfun7. 
y <- c(3, 5.3, 6.4, 9.1, 14.1, 15.4, 18.1, 15.3, 14, 5.8, 7.3, 14.4)
x <- c(1, 1.5, 2,   3,   4,    5,    6,    5,    4,  1,   2,   4.5)
d <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0)
### first we estimate beta, the MLE
lm.wfit(x=cbind(rep(1,12),x), y=y, w=WKM(x=y, d=d)$jump[rank(y)])$coef
## you should get 1.392885 and 2.845658
## then define myfun7 with the MLE value
myfun7 <- function(y, xmat) {
temp1 <- y - ( 1.392885 +  2.845658 * xmat)
return( cbind( temp1, xmat*temp1) )
}
## now test 
el.cen.EM2(y,d, fun=myfun7, mu=c(0,0), xmat=x)
## we should get, Pval = 1 , as the MLE should.
## for other values of (a, b) inside myfun7, you get other Pval
##
rqfun1 <- function(y, xmat, beta, tau = 0.5) {
temp1 <- tau - (1-myfun55(y-beta*xmat))
return(xmat * temp1)
}
myfun55 <- function(x, eps=0.001){
u <- x*sqrt(5)/eps
INDE <- (u < sqrt(5)) & (u > -sqrt(5))
u[u >= sqrt(5)] <- 0
u[u <= -sqrt(5)] <- 1
y <- 0.5 - (u - (u)^3/15)*3/(4*sqrt(5))
u[ INDE ] <- y[ INDE ]
return(u)
}
## myfun55 is a smoothed indicator fn. 
## eps should be between (1/sqrt(n), 1/n^0.75) [Chen and Hall]
el.cen.EM2(x=y,d=d,xc=1:12,fun=rqfun1,mu=0,xmat=x,beta=3.08,tau=0.44769875)
## default tau=0.5 
el.cen.EM2(x=y,d=d,xc=1:12,fun=rqfun1,mu=0,xmat=x,beta=3.0799107404)
###################################################
### next 2 examples are testing the mean/median residual time
###################################################
mygfun <- function(s, age, muage) {as.numeric(s >= age)*(s-(age+muage))}
mygfun2 <- function(s, age, Mdage) 
          {as.numeric(s <= (age+Mdage)) - 0.5*as.numeric(s <= age)}
\dontrun{
library(survival) 
time <- cancer$time
status <- cancer$status-1
###for mean residual time 
el.cen.EM2(x=time, d=status, fun=mygfun, mu=0, age=365.25, muage=234)$Pval
el.cen.EM2(x=time, d=status, fun=mygfun, mu=0, age=365.25, muage=323)$Pval
### for median resudual time
el.cen.EM2(x=time, d=status, fun=mygfun2, mu=0.5, age=365.25, Mdage=184)$Pval
el.cen.EM2(x=time, d=status, fun=mygfun2, mu=0.5, age=365.25, Mdage=321)$Pval
}
\dontrun{
#### For right censor only data (Kaplan-Meier) we can use this function to get a faster computation
#### by calling the kmc 0.2-2 package.
el.cen.R <- function (x, d, xc = 1:length(x), fun, mu, error = 1e-09, ...)
{
xvec <- as.vector(x)
d <- as.vector(d)
mu <- as.vector(mu)
xc <- as.vector(xc)
n <- length(d)
if (length(xvec) != n)
stop("length of d and x must agree")
if (length(xc) != n)
stop("length of xc and d must agree")
if (n <= 2 * length(mu) + 1)
stop("Need more observations")
if (any((d != 0) & (d != 1) ))
stop("d must be 0(right-censored) or 1(uncensored)")
if (!is.numeric(xvec))
stop("x must be numeric")
if (!is.numeric(mu))
stop("mu must be numeric")

funx <- as.matrix(fun(xvec, ...))
pp <- ncol(funx)
if (length(mu) != pp)
stop("length of mu and ncol of fun(x) must agree")
temp <- Wdataclean5(z = xvec, d, zc = xc, xmat = funx)
x <- temp$value
d <- temp$dd
w <- temp$weight
funx <- temp$xxmat
d[length(d)] <- 1
xd1 <- x[d == 1]
if (length(xd1) <= 1)
stop("need more distinct uncensored obs.")
funxd1 <- funx[d == 1, ]
xd0 <- x[d == 0]
wd1 <- w[d == 1]
wd0 <- w[d == 0]
m <- length(xd0)

pnew <- NA
num <- NA
if (m > 0) {
gfun <- function(x) { return( fun(x, ...) - mu ) }
temp <- kmc.solve(x=x, d=d, g=list(gfun))
logel <- temp$loglik.h0
logel00 <- temp$loglik.null
lam <- - temp$lambda
}
if (m == 0) {
num <- 0
temp6 <- el.test.wt2(x = funxd1, wt = wd1, mu)
pnew <- temp6$prob
lam <- temp6$lambda
logel <- sum(wd1 * log(pnew))
logel00 <- sum(wd1 * log(wd1/sum(wd1)))
}
tval <- 2 * (logel00 - logel)
list(loglik = logel, times = xd1, prob = pnew, lam = lam,
iters = num, `-2LLR` = tval, Pval = 1 - pchisq(tval,
df = length(mu)))
}

}
}
\keyword{nonparametric}
\keyword{survival}
\keyword{htest}

