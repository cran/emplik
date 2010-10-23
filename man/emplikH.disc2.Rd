\name{emplikH.disc2}
\alias{emplikH.disc2}
\title{Two sample empirical likelihood ratio for discrete hazards
with right censored, left truncated data, one parameter.}
\usage{
emplikH.disc2(x1, d1, y1= -Inf, x2, d2, y2 = -Inf,
              theta, fun1, fun2, tola = 1e-6, maxi, mini)
}
\description{
Use empirical likelihood ratio and Wilks theorem to test 
the null hypothesis that
\deqn{ 
\int{f_1(t) I_{[dH_1 <1]} \log(1-dH_1(t))} -
\int{f_2(t) I_{[dH_2 <1]} \log(1-dH_2(t))} = \theta 
}
where \eqn{H_*(t)} is the (unknown) discrete cumulative
hazard function; \eqn{f_*(t)} can be any predictable  
functions of \eqn{t}. 
\eqn{\theta} is the parameter. The given value of \eqn{\theta}
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
         \code{fun1(x)} must be able to take a vector input \code{x}.}
    \item{fun2}{ same as fun1.}
    \item{tola}{an optional positive real number, the tolerance of
       iteration error in solve the non-linear equation needed in constrained 
        maximization.}
    \item{theta}{a given real number. for Ho constraint. }
    \item{maxi}{ upper bound for lambda, usually positive. }
    \item{mini}{ lower bound for lambda, usually negative. }
}
\value{
    A list with the following components:
    \item{times}{the location of the hazard jumps.}
    \item{wts}{the jump size of hazard function at those locations.}
    \item{lambda}{the final value of the Lagrange multiplier.}
    \item{"-2LLR"}{The -2Log Likelihood ratio.}
    \item{Pval}{P-value}
    \item{niters}{number of iterations used}
}
\details{
The log likelihood been maximized is the binomial likelihood:
\deqn{ 
\sum D_{1i} \log w_i  + (R_{1i}-D_{1i}) \log [1-w_i] + 
\sum D_{2j} \log v_j  + (R_{2j}-D_{2j}) \log [1-v_j] 
}
where \eqn{w_i = \Delta H_1(t_i)} is the jump 
of the cumulative hazard function at \eqn{t_i}, 
\eqn{D_{1i}} is the number of failures
observed at \eqn{t_i}, \eqn{R_{1i}} is 
the number of subjects at risk at
time \eqn{t_i}.

For discrete distributions, the jump size of the cumulative hazard at
the last jump is always 1. We have to exclude this jump from the 
summation in the constraint calculation
since \eqn{ \log( 1- dH(\cdot))} do not make sense.

The constants \code{theta} must be inside the so called
feasible region for the computation to continue. This is similar to the
requirement that in testing the value of the mean, the value must be
inside the convex hull of the observations.
It is always true that the NPMLE values are feasible. So when the
computation stops, try move the \code{theta} closer
to the NPMLE.  When the computation stops, the -2LLR should have value
infinite.

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
emplikH.disc2(x1=xfemale, d1=dfemale, y1=yfemale, 
 x2=xmale, d2=dmale, y2=ymale, theta=0.2, fun1=fun1, fun2=fun1, maxi=4, mini=-10)
######################################################
### You should get "-2LLR" = 1.511239 and a lot more other outputs.
########################################################
emplikH.disc2(x1=xfemale, d1=dfemale, y1=yfemale, 
 x2=xmale, d2=dmale, y2=ymale, theta=0.25, fun1=fun1, fun2=fun1, maxi=4, mini=-5)
########################################################
### This time you get "-2LLR" = 1.150098 etc. etc.
##############################################################
}
}
\keyword{nonparametric}
\keyword{survival}
\keyword{htest}
