\name{el.trun.test}
\alias{el.trun.test}
\title{Empirical likelihood ratio for mean with left truncated data}
\usage{
el.trun.test(y,x,fun=function(t){t},mu,maxit=20,error=1e-9)
}
\description{
This program uses EM algorithm to compute the maximized
(wrt \eqn{p_i}) empirical
log likelihood function for left truncated data with 
the MEAN constraint:
\deqn{ \sum  p_i f(x_i)  = \int f(t) dF(t) = \mu ~. }
Where \eqn{p_i = \Delta F(x_i)} is a probability.
\eqn{\mu} is a given constant. 
It also returns those \eqn{p_i} and the \eqn{p_i} without
constraint, the Lynden-Bell estimator.

The log likelihood been maximized is
\deqn{ \sum_{i=1}^n \log \frac{\Delta F(x_i)}{1-F(y_i)} .}
}
\arguments{
    \item{y}{a vector containing the left truncation times.}
    \item{x}{a vector containing the survival times. truncation means x>y.}
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
    \item{"-2LLR"}{the maximized empirical log likelihood ratio
                  under the constraint.}
    \item{NPMLE}{jumps of NPMLE of CDF at ordered x.}
    \item{NPMLEmu}{same jumps but for constrained NPMLE.}
}
\details{

This implementation is all in R and have several for-loops in it. 
A faster version would use C to do the for-loop part. But it seems
faster enough and is easier to port to Splus.

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

}
\author{ Mai Zhou }
\references{
    Zhou, M. (2002). 
        Computing censored empirical likelihood ratio 
        by EM algorithm. 
    \emph{Tech Report, Univ. of Kentucky, Dept of Statistics}

    Li, G. (1995)
        Nonparametric likelihood ratio estimation of probabilities
        for truncated data. \emph{JASA} 90, 997-1003.
    
    Turnbull (1976)
       The empirical distribution function with arbitrarily grouped, censored
and truncated data. \emph{JRSS} B 38, 290-295.

}
\examples{
## example with tied observations
vet <- c(30, 384, 4, 54, 13, 123, 97, 153, 59, 117, 16, 151, 22, 56, 21, 18,
         139, 20, 31, 52, 287, 18, 51, 122, 27, 54, 7, 63, 392, 10)
vetstart <- c(0,60,0,0,0,33,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
el.trun.test(vetstart, vet, mu=80, maxit=15)
}
\keyword{nonparametric}
\keyword{survival}
\keyword{htest}

