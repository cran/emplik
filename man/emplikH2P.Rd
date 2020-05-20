\name{emplikH2P}
\alias{emplikH2P}
\title{Return Poisson Empirical likelihood ratio for the given lambda, with 2-sample right censored data}
\usage{
emplikH2P(lambda, x1, d1, x2, d2, fun1, fun2, CIforTheta=FALSE)
}
\description{Compute the Poisson empirical likelihood ratio for the given tilt parameter lambda.
Most useful when construct Wilks confidence intervals. 
The null hypothesis or constraint is defined by the parameter \eqn{\theta}, where 
\deqn{\int fun1(t) dH_1(t) - \int fun2(t) dH_2(t) = \theta }.

If the lambda value set to zero, you get the Nelson-Aalen (NPMLE) and output -2LLR =0.
For other values of lambda, this is not scaled (like in one sample) since there are two samples and it can be
confusing which sample size to use. So here, the lambda is larger (with a scale of sample size) than those in 
one sample: \code{emplikH1P}.
 
Where \eqn{H_1(t)} is the unknown
cumulative hazard
function of sample one; \eqn{H_2(t)} is the cumulative function of sample two; 
\eqn{fun1(t)} can be any given function. 
In the future, the function \eqn{fun1} may replaced by the vector of \eqn{fun(x)}, 
since this is more flexible. Same for \eqn{fun2}.

Input data can be right censored. If no censoring, set \code{d1=rep(1, length(x1))}, and/or \code{d2=rep(1, length(x2))}.
}
\arguments{
\item{lambda}{a scalar. Can be positive or negative. The amount of tiling.}
   \item{x1}{a vector of the censored survival times. sample one.}
   \item{d1}{a vector of the censoring indicators, 1-uncensor; 0-right censor.}
   \item{x2}{a vector of the censored survival times. sample two.}
   \item{d2}{a vector of the censoring indicators, 1-uncensor; 0-right censor.}
   \item{fun1}{a left continuous (weight) function used to calculate 
       the weighted hazard in the parameter \eqn{\theta}. \code{fun1} must be able 
       to take a vector input. See example below.}
    \item{fun2}{Ditto.}
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
EL1 = \sum_{i=1}^n \delta_i \log (dH_1(x_i))  - [ H_1(x_i) ]  ~,
}
(similarly defined for sample 2) and the final EL is the sum of EL1 and EL2.

If there are ties in the data that are resulted from rounding,
you may break the tie by adding a different tiny number to the tied
observation(s). For example: 2, 2, 2, change to 2.00001, 2.00002, 2.00003.
If those are true ties 
(thus the true distribution must be discrete)
we recommend use \code{emplikH2B}. 
}
\value{
    A list with the following components: 
	\item{"-2LLR"}{the -2Log Empirical Likelihood ratio, Poisson version.}
    \item{lambda}{The tilt parameter used. It is also the Lagrange multiplier.}
	\item{"-2LLR(sample1)"}{the -2Log EL ratio, sample 1, Poisson version. "-2LLR" = -2LLR(sample1) + -2LLR(sample2)}
    \item{HazDiff}{Average hazard for the constrained hazard integral, if CIforTheta =TRUE.}
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
