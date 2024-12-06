\name{emplikH1P}
\alias{emplikH1P}
\title{Return Poisson Empirical likelihood ratio for the given lambda, with right censored data}
\usage{
emplikH1P(lambda, x, d, fung, CIforTheta=FALSE)
}
\description{Compute the Poisson empirical likelihood ratio for the given tilt parameter lambda. 
Most useful for the construction of Wilks confidence intervals.
The null hypothesis or constraint is defined by the parameter \eqn{\theta}, where 
\deqn{\int fung(t) dH(t) = \theta }

Where \eqn{H(t)} is the unknown
cumulative hazard
function; \eqn{fung(t)} can be any given function. 

In the future, the function \eqn{fung} may replaced by the vector of \eqn{fung(x)}, 
since this is more flexible.

Input data can be right censored. If no censoring, set \code{d=rep(1, length(x))}. 
}
\arguments{
\item{lambda}{a scalar. Can be positive or negative. The amount of tiling.}
   \item{x}{a vector of the censored survival times.}
   \item{d}{a vector of the censoring indicators, 1-uncensor; 0-right censor.}
   \item{fung}{a left continuous (weight) function used to calculate 
       the weighted hazard in the parameter \eqn{\theta}. \code{fung} must be able 
       to take a vector input. See example below.}
    \item{CIforTheta}{an optional logical value. Default to FALSE. If set to TRUE,
	         will return the integrated hazard value for the given lambda.}
}
\details{
This function is for calculate lambda confidence intervals for \eqn{\theta}.
 
This function is designed for the case where the 
true distribution should be continuous.
So there should be no tie in the data.

The log empirical likelihood used here is the `Poisson' version empirical likelihood:
\deqn{  
\sum_{i=1}^n \delta_i \log (dH(x_i))  - [ H(x_i) ]  ~.
}

If there are ties in the data that are resulted from rounding,
you may want to break the tie by adding a different tiny number to the tied
observation(s). For example: 2, 2, 2, change to 2.00001, 2.00002, 2.00003.
If those are true ties 
(thus the true distribution must be discrete)
we recommend to use \code{emplikH1B} instead. 
}
\value{
    A list with the following components:
    \item{times}{the location of the hazard jumps.}
    \item{wts}{the jump size of hazard function at those locations.}
    \item{lambda}{the Lagrange multiplier.}
    \item{"-2LLR"}{the -2Log Empirical Likelihood ratio, Poisson version.}
    \item{MeanHaz}{The theta defined above, the hazard integral, if CIforTheta =TRUE.}
}
\author{ Mai Zhou } 
\references{
    Pan, X. and Zhou, M. (2002),
    ``Empirical likelihood in terms of hazard for censored data''. 
    \emph{Journal of Multivariate Analysis} \bold{80}, 166-188.
}
\examples{
## fun <- function(x) { as.numeric(x <= 6.5) }
## emplikH1.test( x=c(1,2,3,4,5), d=c(1,1,0,1,1), theta=2, fun=fun) 
## fun2 <- function(x) {exp(-x)}  
## emplikH1.test( x=c(1,2,3,4,5), d=c(1,1,0,1,1), theta=0.2, fun=fun2) 
}
\keyword{nonparametric}
\keyword{survival}
