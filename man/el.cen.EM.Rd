\name{el.cen.EM}
\alias{el.cen.EM}
\title{Empirical likelihood ratio for mean 
with right, left or doubly censored data, by EM algorithm}
\usage{
el.cen.EM(x,d,fun=function(x){x},mu,maxit=25,error=1e-9)
}
\description{
This program uses EM algorithm to compute the maximized 
(wrt \eqn{p_i}) empirical
log likelihood function for right, left or doubly censored data with 
the MEAN constraint:
\deqn{ \sum_{d_i=1}  p_i f(x_i)  = \int f(t) dF(t) = \mu ~. }
Where \eqn{p_i = \Delta F(x_i)} is a probability,
\eqn{d_i} is the censoring indicator, 1(uncensored), 0(right censored),
2(left censored). 
The \eqn{d} for the largest observation is always taken to be 1.
\eqn{\mu} is a given constant. 
}
\arguments{
    \item{x}{a vector containing the observed survival times.}
    \item{d}{a vector containing the censoring indicators, 
           1-uncensored; 0-right censored; 2-left censored.}
    \item{fun}{a continuous (weight) function used to calculate
         the mean as in \eqn{H_0}.
         \code{fun(x)} must be able to take a vector input \code{x}.
         Default to the identity function \eqn{f(x)=x}.}
    \item{mu}{a real number used in the constraint, mean value of f(X).}
    \item{error}{an optional positive real number specifying the tolerance of
       iteration error. This is the bound of the
       \eqn{L_1} norm of the differnence of two successive weights.}
    \item{maxit}{an optional integer, used to control maximum number of
             iterations. }
}
\value{
    A list with the following components:
    \item{loglik}{the maximized empirical log likelihood under the constraint.}
    \item{times}{locations of CDF that have positive mass.}
    \item{prob}{the jump size of CDF function at those locations.}
}
\details{

This implementation is all in R and have several for-loops in it. 
A better version would use C to do the for-loop part.

We only return the log likelihood, not the -2 log likelihood ratio.
You have to plot a curve with many values of the parameter to
find out where is the place the log likelihood becomes maximum.
And from there you can get -2 log likelihood ratio between
the maximum location and your current parameter in Ho.

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
\eqn{p_i^0} taken to be the jumps of NPMLE. 
Or use a different \code{fun}. 

}
\author{ Mai Zhou }
\references{
    Zhou (2002). 
        Computing censored empirical likelihood ratio 
        by EM algorithm. 
    \emph{Tech Report, Univ. of Kentucky, Dept of Statistics}

    Murphy and van der Varrt (1997)
         Semiparametric likelihood ratio inference.
         \emph{Ann. Statist.} \bold{ 25}, 1471-1509.
}
\examples{
## example with tied observations
x <- c(1, 1.5, 2, 3, 4, 5, 6, 5, 4, 1, 2, 4.5)
d <- c(1,   1, 0, 1, 0, 1, 1, 1, 1, 0, 0,   1)
el.cen.EM(x,d,mu=3.5)
}
\keyword{nonparametric}
\keyword{survival}
\keyword{htest}

