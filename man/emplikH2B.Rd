\name{emplikH2B}
\alias{emplikH2B}
\title{Return binomial empirical likelihood ratio for the given lambda, with 2-sample right censored data}
\usage{
emplikH2B(lambda, x1, d1, x2, d2, fun1, fun2, CIforTheta=FALSE)
}
\description{Compute the binomial empirical likelihood ratio for the given tilt parameter lambda. 
Most useful for construct Wilks confidence intervals.
The null hypothesis or constraint is defined by the parameter \eqn{\theta}, where 
\deqn{
\int fun1(t) d \log(1-H_1(t)) - \int fun2(t) d \log (1-H_2(t)) = \theta }.

If the lambda=0, you get the Nelson-Aalen (NPMLE) and output -2LLR =0. Otherwise the lambda is not scaled
(as in one sample case). Since there are two sample sizes. It can be confusing which sample size to use for scale.
So the lambda here is larger than those in one sample by a sclae of (either?) sample size.

 
Where \eqn{H_1(t)} and \eqn{H_2(t)} are the unknown
cumulative hazard function for sample 1/2; \eqn{fun1(t)} and \eqn{fun2(t)} can be any given function. 
It can even be random, just need to be predictable.
In the future, the input function \eqn{fun} may replaced by the vector of \eqn{fun(x)}, 
since this is more flexible.

Input data can be right censored. If no censoring, set \code{d1=rep(1, length(x1))}, and/or \code{d2=rep(1, length(x2))}.
}
\arguments{
\item{lambda}{a scalar. Can be positive or negative. The amount of tiling.}
   \item{x1}{a vector of the censored survival times. sample 1}
   \item{d1}{a vector of the censoring indicators, 1-uncensor; 0-right censor.}
   \item{x2}{a vector of the censored survival times. sample 2}
   \item{d2}{a vector of the censoring indicators, 1-uncensor; 0-right censor.}
   \item{fun1}{a left continuous (weight) function used to calculate 
       the weighted hazard in the parameter \eqn{\theta}. \code{fun1} must be able 
       to take a vector input. See example below.}
    \item{fun2}{Ditto}
    \item{CIforTheta}{an optional logical value. Default to FALSE. If set to TRUE,
	         will return the integrated hazard value for the given lambda.}
}
\details{
This function is used to calculate lambda confidence interval (Wilks type) for \eqn{\theta}.

This function is designed for the case where the 
true distribution should be discrete. Ties in the data are OK.

The log empirical likelihood used here is the `binomial' version empirical likelihood:
\deqn{  
log EL1 = \sum_{i=1}^n \delta_i \log (dH(x_i))  + (R_i - \delta_i)\log [1- dH(x_i) ]  ~,
}
(similarly defined for sample 2) and the overall log EL = log EL1 + log EL2.

}
\value{
    A list with the following components: 
	\item{"-2LLR"}{the -2Log Empirical Likelihood ratio, binomial version.}
	\item{lambda}{the input lambda. The tilt. The Lagrange multiplier.}
    \item{times1}{the location of the hazard jumps. sample 1.} 
	\item{times2}{the location of the hazard jumps. sample 2.}
    \item{wts1}{the jump size of hazard function at those locations.}
    \item{wts2}{the jump size of hazard function at those locations.}
    \item{HazDiff2}{Difference of two hazard integrals. theta defined above.}
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
