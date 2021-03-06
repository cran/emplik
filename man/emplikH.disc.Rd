\name{emplikH.disc}
\alias{emplikH.disc}
\title{Empirical likelihood ratio for discrete hazard
with right censored, left truncated data}
\usage{
emplikH.disc(x, d, y= -Inf, K, fun, tola=.Machine$double.eps^.25, theta)
}
\description{
Use empirical likelihood ratio and Wilks theorem to test 
the null hypothesis that
\deqn{ 
\sum_i[f(x_i, \theta) \log(1- dH(x_i))] = K  
}
where \eqn{H(t)} is the (unknown) discrete cumulative
hazard function; \eqn{f(t,\theta)} can be any predictable
function of \eqn{t}. 
\eqn{\theta} is the parameter of the function 
and \code{K} is a given constant. 
The data can be right censored and left truncated.

When the given constants \eqn{\theta} and/or \code{K} are too far
away from the NPMLE, there will be no hazard function satisfy this 
constraint and the minus 2Log empirical likelihood ratio
will be infinite. In this case the computation will stop.
}
\arguments{
    \item{x}{a vector, the observed survival times.}
    \item{d}{a vector, the censoring indicators, 1-uncensor; 0-censor.}
    \item{y}{optional vector, the left truncation times.}
    \item{K}{a real number used in the constraint, sum to this value.}
    \item{fun}{a left continuous (weight) function used to calculate
         the weighted discrete hazard in \eqn{H_0}. 
         \code{fun(x, theta)} must be able to take a
       vector input \code{x}, and a parameter \code{theta}. }
    \item{tola}{an optional positive real number specifying the tolerance of
       iteration error in solve the non-linear equation needed in constrained 
        maximization.}
    \item{theta}{a given real number used as the parameter of the 
	      function \eqn{f}.}
}
\value{
    A list with the following components:
    \item{times}{the location of the hazard jumps.}
    \item{wts}{the jump size of hazard function at those locations.}
    \item{lambda}{the final value of the Lagrange multiplier.}
    \item{"-2LLR"}{The discrete -2Log Likelihood ratio.}
    \item{Pval}{P-value}
    \item{niters}{number of iterations used}
}
\details{
The log likelihood been maximized is the `binomial' empirical likelihood:
\deqn{ \sum D_i \log w_i  + (R_i-D_i) \log [1-w_i] }
where \eqn{w_i = \Delta H(t_i)} is the jump 
of the cumulative hazard function, \eqn{D_i} is the number of failures
observed at \eqn{t_i}, \eqn{R_i} is the number of subjects at risk at
time \eqn{t_i}.

For discrete distributions, the jump size of the cumulative hazard at
the last jump is always 1. We have to exclude this jump from the 
summation since \eqn{ \log( 1- dH(\cdot))} do not make sense.

The constants \code{theta} and \code{K} must be inside the so called
feasible region for the computation to continue. This is similar to the
requirement that in testing the value of the mean, the value must be
inside the convex hull of the observations.
It is always true that the NPMLE values are feasible. So when the
computation stops, try move the \code{theta} and \code{K} closer
to the NPMLE.  When the computation stops, the -2LLR should have value
infinite.

In case you do not need the \code{theta} in the definition of the
function \eqn{f}, you still need to formally define your \code{fun} function
with a \code{theta} input, just to match the arguments.
}
\author{ Mai Zhou }
\references{
    Fang, H. (2000). Binomial Empirical Likelihood Ratio Method in
	   Survival Analysis. Ph.D. Thesis, Univ. of Kentucky, Dept
 	   of Statistics.

    Zhou and Fang (2001). 
    ``Empirical likelihood ratio for 2 sample problem 
                       for censored data''. 
    \emph{Tech Report, Univ. of Kentucky, Dept of Statistics}

    Zhou, M. and Fang, H. (2006). 
    A comparison of Poisson and binomial empirical likelihood.
\emph{Tech Report, Univ. of Kentucky, Dept of Statistics}
}
\examples{
fun4 <- function(x, theta) { as.numeric(x <= theta) }
x <- c(1, 2, 3, 4, 5, 6, 5, 4, 3, 4, 1, 2.4, 4.5)
d <- c(1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1)
# test if -H(4) = -0.7 
emplikH.disc(x=x,d=d,K=-0.7,fun=fun4,theta=4)
# we should get "-2LLR" 0.1446316  etc....
y <- c(-2,-2, -2, 1.5, -1)
emplikH.disc(x=x,d=d,y=y,K=-0.7,fun=fun4,theta=4)
}
\keyword{nonparametric}
\keyword{survival}
\keyword{htest}
