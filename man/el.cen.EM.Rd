\name{el.cen.EM}
\alias{el.cen.EM}
\title{Empirical likelihood ratio for mean 
with right, left or doubly censored data, by EM algorithm}
\usage{
el.cen.EM(x,d,wt=rep(1,length(d)),fun=function(t){t},mu,maxit=25,error=1e-9,...)
}
\description{
This program uses EM algorithm to compute the maximized 
(wrt \eqn{p_i}) empirical
log likelihood function for right, left or doubly censored data with 
the MEAN constraint:
\deqn{ \sum_{d_i=1}  p_i f(x_i)  = \int f(t) dF(t) = \mu . }
Where \eqn{p_i = \Delta F(x_i)} is a probability,
\eqn{d_i} is the censoring indicator, 1(uncensored), 0(right censored),
2(left censored). 
It also returns those \eqn{p_i}. 

The empirical log likelihood been maximized is
\deqn{ \sum_{d_i=1} \log \Delta F(x_i) + \sum_{d_i=0} \log [1-F(x_i)] 
    + \sum_{d_i=2}  \log F(x_i) . } 
}
\arguments{
    \item{x}{a vector containing the observed survival times.}
    \item{d}{a vector containing the censoring indicators, 
           1-uncensored; 0-right censored; 2-left censored.}
    \item{wt}{a weight vector (case weight). positive. same length as d}
    \item{fun}{a left continuous (weight) function used to calculate
         the mean as in \eqn{H_0}.
         \code{fun(t)} must be able to take a vector input \code{t}.
         Default to the identity function \eqn{f(t)=t}.} 
    \item{mu}{a real number used in the constraint, the mean value of \eqn{f(X)}.}
    \item{maxit}{an optional integer, used to control maximum number of
             iterations. }
    \item{error}{an optional positive real number specifying the tolerance of
       iteration error. This is the bound of the
       \eqn{L_1} norm of the difference of two successive weights.}
   \item{...}{additional arguments, if any, to pass to fun.}
}
\value{
    A list with the following components:
    \item{loglik}{the maximized empirical log likelihood under the constraint.}
    \item{times}{locations of CDF that have positive mass.}
    \item{prob}{the jump size of CDF at those locations.}
    \item{"-2LLR"}{If available, it is Minus two times the 
                   Empirical Log Likelihood Ratio.
                   Should be approximately chi-square distributed under Ho.}
    \item{Pval}{The P-value of the test, using chi-square approximation.}
    \item{lam}{The Lagrange multiplier. Added 5/2007.}
}
\details{

This implementation is all in R and have several for-loops in it. 
A faster version would use C to do the for-loop part.
But this version seems faster enough and is easier to port to Splus.

We return the log likelihood all the time. Sometimes, (for right censored
and no censor case) we also return the -2 log likelihood ratio.
In other cases, you have to plot a curve with many values of the 
parameter, mu, to
find out where is the place the log likelihood becomes maximum.
And from there you can get -2 log likelihood ratio between
the maximum location and your current parameter in Ho.

In order to get a proper distribution as NPMLE, we automatically
change the \eqn{d} for the largest observation to 1
(even if it is right censored), similar for the left censored, 
smallest observation.
\eqn{\mu} is a given constant. 
When the given constants \eqn{\mu} is too far
away from the NPMLE, there will be no distribution
satisfy the constraint.
In this case the computation will stop.
The -2 Log empirical likelihood ratio
should be infinite. 

The constant \code{mu} must be inside 
\eqn{( \min f(x_i) , \max f(x_i) ) }
for the computation to continue. 
It is always true that the NPMLE values are feasible. So when the
computation stops, try move the \code{mu} closer
to the NPMLE --- 
\deqn{ \sum_{d_i=1} p_i^0 f(x_i) } 
\eqn{p_i^0} taken to be the jumps of the NPMLE of CDF. 
Or use a different \code{fun}. 

Difference to the function \code{el.cen.EM2}: here duplicate (input) observations 
are collapsed (with weight 2, 3, ... etc.) but those
will stay separate by default in the \code{el.cen.EM2}. This will lead to
a different \code{loglik} value. But the \code{-2LLR} value should be same
in either version. 

}
\author{ Mai Zhou }
\references{

Zhou, M. (2005). Empirical likelihood ratio with arbitrary censored/truncated data by EM algorithm. 
 \emph{Journal of Computational and Graphical Statistics}, 643-656.

    Murphy, S. and van der Vaart (1997)
         Semiparametric likelihood ratio inference.
         \emph{Ann. Statist.} \bold{ 25}, 1471-1509.
}
\examples{
## example with tied observations
x <- c(1, 1.5, 2, 3, 4, 5, 6, 5, 4, 1, 2, 4.5)
d <- c(1,   1, 0, 1, 0, 1, 1, 1, 1, 0, 0,   1)
el.cen.EM(x,d,mu=3.5)
## we should get "-2LLR" = 1.2466....
myfun5 <- function(x, theta, eps) {
u <- (x-theta)*sqrt(5)/eps 
INDE <- (u < sqrt(5)) & (u > -sqrt(5)) 
u[u >= sqrt(5)] <- 0 
u[u <= -sqrt(5)] <- 1 
y <- 0.5 - (u - (u)^3/15)*3/(4*sqrt(5)) 
u[ INDE ] <- y[ INDE ] 
return(u)
}
el.cen.EM(x, d, fun=myfun5, mu=0.5, theta=3.5, eps=0.1)
## example of using wt in the input. Since the x-vector contain
## two 5 (both d=1), and two 2(both d=0), we can also do
xx <- c(1, 1.5, 2, 3, 4, 5, 6, 4, 1, 4.5)
dd <- c(1,   1, 0, 1, 0, 1, 1, 1, 0,   1)
wt <- c(1,   1, 2, 1, 1, 2, 1, 1, 1,   1)
el.cen.EM(x=xx, d=dd, wt=wt, mu=3.5)
## this should be the same as the first example.
}
\keyword{nonparametric}
\keyword{survival}
\keyword{htest}

