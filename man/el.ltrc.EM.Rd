\name{el.ltrc.EM}
\alias{el.ltrc.EM}
\title{Empirical likelihood ratio for mean 
with left truncated and right censored data, by EM algorithm}
\usage{
el.ltrc.EM(y,x,d,fun=function(t){t},mu,maxit=30,error=1e-9)
}
\description{
This program uses EM algorithm to compute the maximized 
(wrt \eqn{p_i}) empirical
log likelihood function for left truncated and right censored data with 
the MEAN constraint:
\deqn{ \sum_{d_i=1}  p_i f(x_i)  = \int f(t) dF(t) = \mu ~. }
Where \eqn{p_i = \Delta F(x_i)} is a probability,
\eqn{d_i} is the censoring indicator, 1(uncensored), 0(right censored). 
The \eqn{d} for the largest observation \eqn{x}, is always (automatically)
changed to 1.  \eqn{\mu} is a given constant. 
This function also returns those \eqn{p_i}. 

The log empirical likelihood function been maximized is
\deqn{\sum_{d_i=1} \log \frac{ \Delta F(x_i)}{1-F(y_i)}  + 
    \sum_{d_i=0} \log \frac{1-F(x_i)}{1-F(y_i)}.}
}
\arguments{
    \item{y}{an optional vector containing the observed truncation times.}
    \item{x}{a vector containing the observed survival times.}
    \item{d}{a vector containing the censoring indicators, 
           1-uncensored; 0-right censored.}
    \item{fun}{a continuous (weight) function used to calculate
         the mean as in \eqn{H_0}.
         \code{fun(t)} must be able to take a vector input \code{t}.
         Default to the identity function \eqn{f(t)=t}.}
    \item{mu}{a real number used in the constraint, mean value of \eqn{f(X)}.}
    \item{error}{an optional positive real number specifying the tolerance of
       iteration error. This is the bound of the
       \eqn{L_1} norm of the difference of two successive weights.}
    \item{maxit}{an optional integer, used to control maximum number of
             iterations. }
}
\value{
    A list with the following components:
    \item{times}{locations of CDF that have positive mass.}
    \item{prob}{the probability of the constrained NPMLE of 
                CDF at those locations.}
    \item{"-2LLR"}{It is Minus two times the 
                   Empirical Log Likelihood Ratio.
                   Should be approximate chi-square distributed under Ho.}
}
\details{

We return the -2 log likelihood ratio, and the constrained
NPMLE of CDF.
The un-constrained NPMLE should be WJT or Lynden-Bell estimator.

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

This implementation is all in R and have several for-loops in it. 
A faster version would use C to do the for-loop part.
(but this version is easier to port to Splus, and seems faster enough).
}
\author{ Mai Zhou }
\references{
    Zhou, M. (2002). 
        Computing censored and truncated empirical likelihood ratio 
        by EM algorithm. 
    \emph{Tech Report, Univ. of Kentucky, Dept of Statistics}

    Turnbbull, B. (1976). The empirical distribution function with
   arbitrarily grouped, censored and truncated data. JRSS B 290-295.
}
\examples{
## example with tied observations
y <- c(0, 0, 0.5, 0, 1, 2, 2, 0, 0, 0, 0, 0 )
x <- c(1, 1.5, 2, 3, 4, 5, 6, 5, 4, 1, 2, 4.5)
d <- c(1,   1, 0, 1, 0, 1, 1, 1, 1, 0, 0,   1)
el.ltrc.EM(y,x,d,mu=3.5)
ypsy <- c(51, 58, 55, 28, 25, 48, 47, 25, 31, 30, 33, 43, 45, 35, 36)
xpsy <- c(52, 59, 57, 50, 57, 59, 61, 61, 62, 67, 68, 69, 69, 65, 76)
dpsy <- c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1 )
el.ltrc.EM(ypsy,xpsy,dpsy,mu=64)
}
\keyword{nonparametric}
\keyword{survival}
\keyword{htest}

