\name{emplikH2.test}
\alias{emplikH2.test}
\title{Empirical likelihood for weighted hazard with 
        right censored, left truncated data}
\usage{
emplikH2.test(x, d, y= -Inf, K, fun, tola=.Machine$double.eps^.25,...)
}
\description{ Use empirical likelihood ratio and Wilks theorem to test 
the null hypothesis that 
\deqn{\int f(t, ... ) dH(t) = K }
with right censored data, where \eqn{ H(t) } is the (unknown)
cumulative hazard function;
\eqn{ f(t, ... )} can be any given left continuous function in \eqn{t};
(of course the integral must be finite). 
}
\arguments{
    \item{x}{a vector containing the observed survival times.}
    \item{d}{a vector of the censoring indicators, 1-uncensor; 0-censor.}
    \item{y}{a vector containing the left truncation times.}
    \item{K}{a real number used in the constraint, i.e. to set the
               weighted integral of hazard to this value.}
    \item{fun}{a left continuous (in \code{t}) weight function
              used to calculate
         the weighted hazard in \eqn{H_0}. \code{fun(t, ... )} must be
         able to take a vector input \code{t}.}
    \item{tola}{an optional positive real number specifying the tolerance of
       iteration error in solve the non-linear equation needed in constrained
        maximization.}
    \item{...}{additional parameter(s) passing along to \code{fun}. This
           allows an implicit function of \code{fun}.  }
}
\value{
    A list with the following components:
    \item{times}{the location of the hazard jumps.}
    \item{wts}{the jump size of hazard function at those locations.}
    \item{lambda}{the Lagrange multiplier.}
    \item{"-2LLR"}{the -2Log Likelihood ratio.}
    \item{Pval}{P-value}
    \item{niters}{number of iterations used}
}
\author{ Mai Zhou } 
\references{
    Pan and Zhou (2002),
    ``Empirical likelihood in terms of cumulative
      hazard for censored data''. 
    \emph{Journal of Multivariate Analysis} \bold{80}, 166-188.
}
\details{
This version works for implicit function \code{ f(t, ...)}.

This function is designed for continuous distributions.
Thus we do not expect tie in the observations. If you believe
the true underlying distribution is continuous but the
sample observations have tie due to rounding, then you might want
to add a small number to the observations to break tie.

The likelihood used here is the `Poisson' version of the likelihood
\deqn{
\prod ( dH(t_i) )^{\delta_i} \exp [-H(t_i)] . 
}

For discrete distributions we recommend
use \code{emplikdisc.test()}.

The constants \code{theta} and \code{K} must be inside the so called
feasible region for the computation to continue. This is similar to the
requirement that in testing the value of the mean, the value must be
inside the convex hull of the observations.
It is always true that the NPMLE values are feasible. So when the
computation stops, that means there is no hazard function satisfy 
the constraint. You may try to move the \code{theta} and \code{K} closer
to the NPMLE.  When the computation stops, the -2LLR should have value
infinite.
}
\examples{
z1<-c(1,2,3,4,5)
d1<-c(1,1,0,1,1)
fun4 <- function(x, theta) { as.numeric(x <= theta) }
emplikH2.test(x=z1,d=d1, K=0.5, fun=fun4, theta=3.5)
#Next, test if H(3.5) = log(2) .
emplikH2.test(x=z1,d=d1, K=log(2), fun=fun4, theta=3.5)
#Next, try one sample log rank test
indi <- function(x,y){ as.numeric(x >= y) }
fun3 <- function(t,z){rowsum(outer(z,t,FUN="indi"),group=rep(1,length(z)))} 
emplikH2.test(x=z1, d=d1, K=sum(0.25* z1), fun=fun3, z=z1) 
}
\keyword{nonparametric}
\keyword{survival}
\keyword{htest}
