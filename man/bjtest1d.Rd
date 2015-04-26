\name{bjtest1d}
\alias{bjtest1d}
\title{Test the Buckley-James estimator by Empirical Likelihood, 1-dim only}
\usage{
bjtest1d(y, d, x, beta)
}
\arguments{
    \item{y}{a vector of length N, containing the censored responses.}
    \item{d}{a vector of either 1's or 0's. d=1 means y is uncensored.
	d=0 means y is right censored. }
    \item{x}{a vector of length N, covariate.}
    \item{beta}{a number. the regression coefficient to 
          be tested in the model y = x beta + epsilon }
}
\description{
    Use the empirical likelihood ratio and Wilks theorem to test if the
    regression coefficient is equal to beta. For 1-dim beta only.

The log empirical likelihood been maximized is
\deqn{ \sum_{d=1} \log \Delta F(e_i) + \sum_{d=0} \log [1-F(e_i)] .}
}
\details{

In the above likelihood, \eqn{ e_i = y_i - x * beta } is the residuals.
 
Similar to \code{bjtest( )}, but only for 1-dim beta.

}
\value{
    A list with the following components:
    \item{"-2LLR"}{the -2 loglikelihood ratio; have approximate chi square 
                  distribution under \eqn{H_o}.}
    \item{logel2}{the log empirical likelihood, under estimating equation.}
    \item{logel}{the log empirical likelihood of the Kaplan-Meier of e's.}
    \item{prob}{the probabilities that max the empirical likelihood 
               under estimating equation constraint.}
}
\references{
  Buckley, J. and James, I. (1979). Linear regression with censored data.
     Biometrika,  \bold{66} 429-36.

  Owen, A. (1990). Empirical likelihood ratio confidence regions. 
    Ann. Statist. \bold{18} 90-120.

  Zhou, M. Li, G. (2004). Empirical likelihood analysis of the
       Buckley-James estimator. Tech. Report.
}
\author{Mai Zhou.}
\examples{
xx <- c(28,-44,29,30,26,27,22,23,33,16,24,29,24,40,21,31,34,-2,25,19)
}
\keyword{nonparametric}
\keyword{htest}
