\name{el.test}
\alias{el.test}
\title{Empirical likelihood ratio test for the means, uncensored data}
\usage{
el.test(x, mu, lam, maxit=25, gradtol=1e-7, 
                 svdtol = 1e-9, itertrace=FALSE)
}
\arguments{
    \item{x}{a matrix or vector containing the data, one row
 	 per observation.}
    \item{mu}{a numeric vector (of length \code{ = ncol(x)})
           to be tested as the mean vector of \code{x}
		above, as \eqn{H_0}.}
    \item{lam}{an optional vector of length \code{ = length(mu)},
            the starting value of Lagrange
	 multipliers, will use \eqn{0} if missing.}
    \item{maxit}{an optional integer to control iteration when solve
	 constrained maximization.}
    \item{gradtol}{an optional real value for convergence test.}
    \item{svdtol}{an optional real value to detect singularity while
		solve equations.}
    \item{itertrace}{a logical value. If the iteration history 
	needs to be printed out.}
}
\description{
    Compute the empirical likelihood ratio with the
	mean vector fixed at mu.
}
\details{
If \code{mu} is in the interior of the convex hull of the
observations \code{x}, then \code{wts} should sum to \code{n}. 
If \code{mu} is outside
the convex hull then \code{wts} should sum to nearly zero, and 
\code{-2LLR}
will be a large positive number.  It should be infinity,
but for inferential purposes a very large number is
essentially equivalent.  If mu is on the boundary of the convex
hull then \code{wts} should sum to nearly k where k is the number of
observations within that face of the convex hull which contains mu.

When \code{mu} is interior to the convex hull, it is typical for
the algorithm to converge quadratically to the solution, perhaps
after a few iterations of searching to get near the solution.
When \code{mu} is outside or near the boundary of the convex hull, then
the solution involves a \code{lambda} of infinite norm.  The algorithm
tends to nearly double \code{lambda} at each iteration and the gradient
size then decreases roughly by half at each iteration.

The goal in writing the algorithm was to have it ``fail gracefully"
when \code{mu} is not inside the convex hull.  The user can 
either leave \code{-2LLR} ``large and positive" or can replace
it by infinity when the weights do not sum to nearly n.
}
\value{
    A list with the following components:
    \item{-2LLR}{the -2 loglikelihood ratio; approximate chisq distribution
                 under \eqn{H_o}.}
    \item{Pval}{the observed P-value by chi-square approximation.}
    \item{lambda}{the final value of Lagrange multiplier.}
    \item{grad}{the gradient at the maximum.}
    \item{hess}{the Hessian matrix.}
    \item{wts}{weights on the observations}
    \item{nits}{number of iteration performed}
}
\references{
    Owen, A. (1990). Empirical likelihood ratio confidence regions. 
    Ann. Statist. \bold{18} 90-120.
}
\author{ Original Splus code by Owen. Adapted to R by Mai Zhou. }
\examples{
x <- matrix(c(rnorm(50,mean=1), rnorm(50,mean=2)), ncol=2,nrow=50)
el.test(x, mu=c(1,2))
## Suppose now we wish to test Ho: 2mu(1)-mu(2)=0, then
y <- 2*x[,1]-x[,2]
el.test(y, mu=0)
xx <- c(28,-44,29,30,26,27,22,23,33,16,24,29,24,40,21,31,34,-2,25,19)
el.test(xx, mu=15)  #### -2LLR = 1.805702
}
\keyword{nonparametric}
\keyword{htest}
