\name{el.cen.kmc1d}
\alias{el.cen.kmc1d}
\title{Empirical likelihood ratio for 1 mean constraint 
with right censored data}
\usage{
el.cen.kmc1d(x, d, fun, mu, tol = .Machine$double.eps^0.5, step=0.001, ...)
}
\description{
This program uses a fast recursive formula to compute the maximized 
(wrt \eqn{p_i}) empirical
log likelihood ratio for right censored data with 
one MEAN constraint:
\deqn{ \sum_{d_i=1}  p_i f(x_i)  = \int f(t) dF(t) = \mu . }
Where \eqn{p_i = \Delta F(x_i)} is a probability,
\eqn{d_i} is the censoring indicator, 1(uncensored), 0(right censored). 
It also returns those \eqn{p_i}. 

The empirical log likelihood been maximized is
\deqn{ \sum_{d_i=1} \log \Delta F(x_i) + \sum_{d_i=0} \log [1-F(x_i)] . } 
}
\arguments{
    \item{x}{a vector containing the observed survival times.}
    \item{d}{a vector containing the censoring indicators, 
           1-uncensored; 0-right censored.}
    \item{fun}{a left continuous (weight) function used to calculate
         the mean as in \eqn{H_0}.
         \code{fun(t)} must be able to take a vector input \code{t}.
         } 
    \item{mu}{a real number used in the constraint, the mean value of \eqn{f(X)}.}
	\item{tol}{a small positive number, for the uniroot error tol.}
	\item{step}{a small positive number, for use in the uniroot function (as interval) to find lambda root. Sometimes 
	            uniroot will find the wrong root or no root, resulting a negative "-2LLR" or NA. Change the step to a different 
				value often can fix this (but not always). Another sign of wrong root is that
				the sum of probabilities not sum to one, or has negative probability values.}
   \item{...}{additional arguments, if any, to pass to fun.}
}
\value{
    A list with the following components:
    \item{loglik}{the maximized empirical log likelihood under the constraint. Note, here the tied observations are
	              not collapsed into one obs. with weight 2 (as in el.cen.EM), so the value may differ from those
				  that do collapse the tied obs. In any case, the -2LLR should not differ (whether collaps or not).}
    \item{times}{locations of CDF that have positive mass.}
    \item{prob}{the jump size of CDF at those locations.}
    \item{"-2LLR"}{If available, it is minus two times the 
                   empirical Log Likelihood Ratio.
                   Should be approximately chi-square distributed under Ho. If you got NA or negative value, then something is wrong,
				   most likely the uniroot has found the wrong root. Suggest: use el.cen.EM2() which uses EM algorithm. It is more stable but slower. }
    \item{Pval}{The P-value of the test, using chi-square approximation.}
    \item{lam}{The Lagrange multiplier.}
}
\details{
This function is similar to the function in package \code{kmc}, but much simpler, i.e.
all implemented in R and only for one mean.
This implementation have two for-loops in R. 
A faster version would use C to do the for-loop part.
But this version seems fast enough and is easier to port to Splus.

We return the log likelihood all the time. Sometimes, (for right censored case) 
we also return the -2 log likelihood ratio.
In other cases, you have to plot a curve with many values of the 
parameter, mu, to
find out where is the place the log likelihood becomes maximum.
And from there you can get -2 log likelihood ratio between
the maximum location and your current parameter in Ho.

The input \code{step} is used in uniroot function to find a root of lambda. Sometimes
a step value may lead to no root or result in a wrong root. You may try several values 
for the step to see. If the probabilities returned do not sum to one, then the lambda root is a wrong root.
We want the root closest to zero.

In order to get a proper distribution as NPMLE, we automatically
change the \eqn{d} for the largest observation to 1
(even if it is right censored).
\eqn{\mu} is a given constant. 
When the given constants \eqn{\mu} is too far
away from the NPMLE, there will be no distribution
satisfy the constraint.
In this case the computation will stop or return something ridiculas, (as negative -2LLR).
The -2 Log empirical likelihood ratio
may be +infinite. 

The constant \code{mu} must be inside 
\eqn{( \min f(x_i) , \max f(x_i) ) } (with uncensored \eqn{x_i})
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
Zhou, M. and Yang, Y. (2015). A recursive formula for the Kaplan-Meier estimator with mean constraints 
and its application to empirical likelihood.
\emph{Computational Statistics}  Vol. 30, Issue 4 pp. 1097-1109.

Zhou, M. (2005). Empirical likelihood ratio with arbitrary censored/truncated data by EM algorithm. 
 \emph{Journal of Computational and Graphical Statistics}, 14(3), 643-656.
}
\examples{
x <- c(1, 1.5, 2, 3, 4.2, 5, 6.1, 5.3, 4.5, 0.9, 2.1, 4.3)
d <- c(1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1)
ff <- function(x) {
    x - 3.7
}
el.cen.kmc1d(x=x, d=d, fun=ff, mu=0)
#######################################
## example with tied observations
x <- c(1, 1.5, 2, 3, 4, 5, 6, 5, 4, 1, 2, 4.5)
d <- c(1,   1, 0, 1, 0, 1, 1, 1, 1, 0, 0,   1)
el.cen.EM(x,d,mu=3.5)
## we should get "-2LLR" = 1.2466....
myfun5 <- function(x, theta, eps) {
u <- (x-theta)*sqrt(5)/eps 
INDE <- (u < sqrt(5)) & (u > -sqrt(5)) 
u[u >= sqrt(5)] <- 0 
u[u <= -sqrt(5)] <- 1 
y <- 0.5 - (u - (u)^3/15)*3/(4*sqrt(5)) 
u[ INDE ] <- y[ INDE ] 
return(u)
}
el.cen.EM(x, d, fun=myfun5, mu=0.5, theta=3.5, eps=0.1)
## example of using wt in the input. Since the x-vector contain
## two 5 (both d=1), and two 2(both d=0), we can also do
xx <- c(1, 1.5, 2, 3, 4, 5, 6, 4, 1, 4.5)
dd <- c(1,   1, 0, 1, 0, 1, 1, 1, 0,   1)
wt <- c(1,   1, 2, 1, 1, 2, 1, 1, 1,   1)
el.cen.EM(x=xx, d=dd, wt=wt, mu=3.5)
## this should be the same as the first example.
}
\keyword{nonparametric}
\keyword{survival}
\keyword{htest}

