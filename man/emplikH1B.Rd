\name{emplikH1B}
\alias{emplikH1B}
\title{Return binomial empirical likelihood ratio for the given lambda, with right censored data}
\usage{
emplikH1B(lambda, x, d, fung, CIforTheta=FALSE)
}
\description{Compute the binomial empirical likelihood ratio for the given tilt parameter lambda. 
Most useful for construct Wilks confidence intervals.
The null hypothesis or constraint is defined by the parameter \eqn{\theta}, where 
\deqn{\int fung(t) d \log (1-H(t)) = \theta }

 
Where \eqn{H(t)} is the unknown
cumulative hazard function; \eqn{fung(t)} can be any given function. 
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
This function is used to calculate lambda confidence interval (Wilks type) for \eqn{\theta}.

This function is designed for the case where the 
true distribution should be discrete. Ties in the data are OK.

The log empirical likelihood used here is the `binomial' version empirical likelihood:
\deqn{  
\sum_{i=1}^n \delta_i \log (dH(x_i))  + (R_i - \delta_i)\log [1- dH(x_i) ]  .
}

}
\value{
    A list with the following components:
    \item{times}{the location of the hazard jumps.}
    \item{jumps}{the jump size of hazard function at those locations.}
    \item{lambda}{the input lambda.}
    \item{"-2LLR"}{the -2Log Likelihood ratio.}
    \item{IntHaz}{The theta defined above, for the given lambda.}
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
