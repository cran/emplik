\name{emplikHs.disc2}
\alias{emplikHs.disc2}
\title{Two sample empirical likelihood ratio for discrete hazards
with right censored, left truncated data. Many constraints.}
\usage{
emplikHs.disc2(x1, d1, y1= -Inf, x2, d2, y2 = -Inf,
          theta, fun1, fun2, maxit=25,tola = 1e-6, itertrace =FALSE)
}
\description{
Use empirical likelihood ratio and Wilks theorem to test 
the null hypothesis that
\deqn{ 
\int{f_1(t) I_{[dH_1 <1]} \log(1-dH_1(t))} -
\int{f_2(t) I_{[dH_2 <1]} \log(1-dH_2(t))} = \theta 
}
where \eqn{H_*(t)} are the (unknown) discrete cumulative
hazard functions; \eqn{f_*(t)} can be any predictable  
functions of \eqn{t}. 
\eqn{\theta} is a vector of parameters (dim=q >= 1). 
The given value of \eqn{\theta}
in these computation are the value to be tested.
The data can be right censored and left truncated.

When the given constants \eqn{\theta} is too far
away from the NPMLE, there will be no hazard function satisfy this 
constraint and the -2 Log empirical likelihood ratio
will be infinite. In this case the computation will stop.
}
\arguments{
    \item{x1}{a vector, the observed survival times, sample 1.}
    \item{d1}{a vector, the censoring indicators, 1-uncensor; 0-censor.}
    \item{y1}{optional vector, the left truncation times.}
    \item{x2}{a vector, the observed survival times, sample 2.}
    \item{d2}{a vector, the censoring indicators, 1-uncensor; 0-censor.}
    \item{y2}{optional vector, the left truncation times.}
    \item{fun1}{a predictable function used to calculate
         the weighted discrete hazard in \eqn{H_0}. 
         \code{fun1(x)} must be able to take a vector input (length n)
         \code{x}, and output a matrix of n x q.}
    \item{fun2}{ Ditto.}
    \item{tola}{an optional positive real number, the tolerance of
       iteration error in solve the non-linear equation needed in constrained 
        maximization.}
    \item{theta}{a given vector of length q. for Ho constraint. }
    \item{maxit}{integer, maximum number of iteration. }
    \item{itertrace}{ Logical. }
}
\value{
    A list with the following components:
    \item{times1}{the location of the hazard jumps in sample 1.}
    \item{times2}{the location of the hazard jumps in sample 2.}
    \item{lambda}{the final value of the Lagrange multiplier.}
    \item{"-2LLR"}{The -2Log Likelihood ratio.}
    \item{"-2LLR(sample1)"}{The -2Log Likelihood ratio for sample 1 only.}
    \item{niters}{number of iterations used}
}
\details{
The log empirical likelihood been maximized is the `binomial empirical likelihood':
\deqn{ 
\sum D_{1i} \log w_i  + (R_{1i}-D_{1i}) \log [1-w_i] + 
\sum D_{2j} \log v_j  + (R_{2j}-D_{2j}) \log [1-v_j] 
}
where \eqn{w_i = \Delta H_1(t_i)} is the jump 
of the cumulative hazard function at \eqn{t_i}, 
\eqn{D_{1i}} is the number of failures
observed at \eqn{t_i}, and \eqn{R_{1i}} is 
the number of subjects at risk at
time \eqn{t_i} (for sample one). Similar for sample two.

For discrete distributions, the jump size of the cumulative hazard at
the last jump is always 1. We have to exclude this jump from the 
summation in the constraint calculation
since \eqn{ \log( 1- dH(\cdot))} do not make sense.
In the likelihood, this term contribute a zero (0*Inf).

This function can handle multiple constraints. So dim( \code{theta}) = q.
The constants \code{theta} must be inside the so called
feasible region for the computation to continue. This is similar to the
requirement that in testing the value of the mean, the value must be
inside the convex hull of the observations.
It is always true that the NPMLE values are feasible. So when the
computation stops, try move the \code{theta} closer
to the NPMLE.  When the computation stops, the -2LLR should have value
infinite.

This code can also be used to compute one sample problems.
You need to artificially supply data for sample two
(with minimal sample size (2q+2)), and supply a function 
\code{fun2} that ALWAYS returns zero (zero vector or zero matrix).
In the output, read the -2LLR(sample1).
}
\author{ Mai Zhou }
\references{
    Zhou and Fang (2001). 
    ``Empirical likelihood ratio for 2 sample problems for censored data''. 
    \emph{Tech Report, Univ. of Kentucky, Dept of Statistics}
}
\examples{
if(require("boot", quietly = TRUE)) {
####library(boot)
data(channing)
ymale <- channing[1:97,2]
dmale <- channing[1:97,5]
xmale <- channing[1:97,3]
yfemale <- channing[98:462,2]
dfemale <- channing[98:462,5]
xfemale <- channing[98:462,3]
fun1 <- function(x) { as.numeric(x <= 960) }
########################################################
emplikHs.disc2(x1=xfemale, d1=dfemale, y1=yfemale, 
 x2=xmale, d2=dmale, y2=ymale, theta=0.25, fun1=fun1, fun2=fun1)
########################################################
### This time you get "-2LLR" = 1.150098 etc. etc.
##############################################################
fun2 <- function(x){ cbind(as.numeric(x <= 960), as.numeric(x <= 860))}
############ fun2 has matrix output ###############
emplikHs.disc2(x1=xfemale, d1=dfemale, y1=yfemale, 
 x2=xmale, d2=dmale, y2=ymale, theta=c(0.25,0), fun1=fun2, fun2=fun2)
################# you get "-2LLR" = 1.554386, etc ###########
}
}
\keyword{nonparametric}
\keyword{survival}
\keyword{htest}
