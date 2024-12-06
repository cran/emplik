\name{BJnoint}
\alias{BJnoint}
\title{The Buckley-James censored regression estimator}
\usage{
BJnoint(x, y, delta, beta0 = NA, maxiter=30, error = 0.00001)
}
\arguments{
    \item{x}{a matrix or vector containing the covariate, one row
 	 per observation.}
    \item{y}{a numeric vector of length N, censored responses. }
    \item{delta}{a vector of length N, delta=0/1 for censored/uncensored.}
    \item{beta0}{an optional vector for starting value of iteration.}
    \item{maxiter}{an optional integer to control iterations.}
    \item{error}{an optional positive value to control iterations.}
}
\description{
    Compute the Buckley-James estimator in the regression model 
\deqn{ y_i = \beta x_i + \epsilon_i } 
with right censored \eqn{y_i}. Iteration method.

}
\details{
This function compute the Buckley-James estimator 
when your model do not have an intercept term.
Of course, if you include a column of 1's in the \code{x} matrix, 
it is also OK with this function and it
is equivalent to having an intercept term.
If your model do have an intercept term, then you could also use the function
\code{bj( )} in the \code{Design} library. It should be more refined 
than \code{BJnoint} in the stopping rule for the iterations. However, the variance
estimator \code{bj( )} provided is not consistent.

This function is included here mainly to produce the estimator value
that may provide some useful information with the function \code{bjtest( )}.
For example you may want to test a beta value near the
Buckley-James estimator. 

}
\value{
    A list with the following components:
    \item{beta}{the Buckley-James estimator.}
    \item{iteration}{number of iterations performed.}
}
\references{
    Buckley, J. and James, I. (1979).  Linear regression with censored data.
   \emph{Biometrika}, \bold{66} 429-36.
   
    Zhou, M. (2016). Empirical Likelihood Method in Survival Analysis. CRC Publishing.
}
\author{ Mai Zhou. }
\examples{
x <- matrix(c(rnorm(50,mean=1), rnorm(50,mean=2)), ncol=2,nrow=50)
## Suppose now we wish to test Ho: 2mu(1)-mu(2)=0, then
y <- 2*x[,1]-x[,2]
xx <- c(28,-44,29,30,26,27,22,23,33,16,24,29,24,40,21,31,34,-2,25,19)
}
\keyword{nonparametric}
\keyword{htest}
