\name{emplikH1.test}
\alias{emplikH1.test}
\title{Empirical likelihood for hazard with right censored, 
                 left truncated data}
\usage{
emplikH1.test(x, d, y= -Inf, theta, fun, tola=.Machine$double.eps^.25)
}
\description{Use empirical likelihood ratio and Wilks theorem to test 
the null hypothesis that 
\deqn{\int f(t) dH(t) = \theta }
with right censored data. Where \eqn{H(t)} is the unknown
cumulative hazard
function; \eqn{f(t)} can be any given left continuous function and
\eqn{\theta} a given constant. In fact, \eqn{f(t)} can even be data
dependent, just have to be `predictable'.
}
\arguments{
   \item{x}{a vector of the observed survival times.}
   \item{d}{a vector of the censoring indicators, 1-uncensor; 0-censor.}
   \item{y}{a vector of the observed left truncation times.}
   \item{theta}{a real number used in the \eqn{H_0} to set the hazard 
               to this value.}
   \item{fun}{a left continuous (weight) function used to calculate 
       the weighted hazard in \eqn{ H_0}. \code{fun} must be able 
       to take a vector input. See example below.}
    \item{tola}{an optional positive real number specifying the tolerance of
       iteration error in solve the non-linear equation needed in constrained
        maximization.}
}
\details{
This function is designed for the case where the 
true distributions are all continuous.
So there should be no tie in the data.

The log empirical likelihood used here is the `Poisson' version likelihood:
\deqn{  
\sum_{i=1}^n \delta_i \log (dH(x_i))  - [ H(x_i) - H(y_i) ]  ~.
}

If there are ties in the data that are resulted from rounding,
you may break the tie by adding a tiny number to the tied
observation(s). If those are true ties 
(thus the true distribution is discrete)
we recommend use \code{emplikdisc.test()}. 

The constant \code{theta} must be inside the so called
feasible region for the computation to continue. This is similar to the
requirement that in testing the value of the mean, the value must be
inside the convex hull of the observations.
It is always true that the NPMLE values are feasible. So when the
computation stops, that means there is no hazard function satisfy
the constraint. You may try to move the \code{theta} value closer
to the NPMLE.  When the computation stops, the -2LLR should have value
infinite.
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
    Pan, X. and Zhou, M. (2002),
    ``Empirical likelihood in terms of hazard for censored data''. 
    \emph{Journal of Multivariate Analysis} \bold{80}, 166-188.
}
\examples{
fun <- function(x) { as.numeric(x <= 6.5) }
emplikH1.test( x=c(1,2,3,4,5), d=c(1,1,0,1,1), theta=2, fun=fun) 
fun2 <- function(x) {exp(-x)}  
emplikH1.test( x=c(1,2,3,4,5), d=c(1,1,0,1,1), theta=0.2, fun=fun2) 
}
\keyword{nonparametric}
\keyword{survival}
\keyword{htest}
