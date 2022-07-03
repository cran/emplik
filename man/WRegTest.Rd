\name{WRegTest}
\alias{WRegTest}
\title{Test the case weighted regression estimator by Empirical Likelihood}
\usage{
WRegTest(x, y, delta, beta0, psifun=function(t){t})
}
\arguments{
    \item{x}{a matrix of size N by q. Random design matrix. }
    \item{y}{a vector of length N, containing the censored responses.}
    \item{delta}{a vector (length N) of either 1's or 0's. 
              delta=1 means y is uncensored;
	      delta=0 means y is right censored. }
    \item{beta0}{a vector of length q. The value of the regression 
           coefficient to be tested in the linear model }.
    \item{psifun}{the estimating function. The definition of it determines
                  the type of estimator under testing. }
}
\description{
    Use the empirical likelihood ratio and Wilks theorem to test if the
    regression coefficient is equal to \code{beta0}, 
       by the case weighted estimation method.

The log empirical likelihood been maximized is
\deqn{ \sum_{d=1} \log \Delta F(y_i) + \sum_{d=0} \log [1-F(y_i)].}
}
\details{

The above likelihood should be understood as the likelihood of the 
censored responses \code{y} and \code{delta}. 

This version can handle the model where beta is a vector (of length q).

The estimation equations used when maximize the 
empirical likelihood is 
\deqn{ 0 = \sum \delta_i \Delta F(Y_i) X_i \psi( Y_i - X_i \beta0 ) }
which was described in detail in the reference below.

For median regression (Least Absolute Deviation) estimator, you should
define the 
\code{psifun} as \eqn{+1, -1} or \eqn{0} when \eqn{t} is \eqn{>0, <0 }
or \eqn{ =0}.

For ordinary least squares estimator, \code{psifun} should be the identity function psifun <- function(t){t}.
}
\value{
    A list with the following components:
    \item{"-2LLR"}{the -2 log likelihood ratio; have approximate chisq 
                  distribution under \eqn{H_0}.}
    \item{P-val}{the p-value using the chi-square approximation.}
}
\references{
    Zhou, M.; Kim, M. and Bathke, A. (2012). 
       Empirical likelihood analysis of the case weighted estimator in 
       heteroscastic AFT model.  
    \emph{Statistica Sinica}, \bold{22}, 295-316.
}
\author{ Mai Zhou. }
\examples{
xx <- c(28,-44,29,30,26,27,22,23,33,16,24,29,24,40,21,31,34,-2,25,19)
}
\keyword{nonparametric}
\keyword{htest}
