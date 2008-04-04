\name{ROCnp}
\alias{ROCnp}
\title{Test the ROC curve by Empirical Likelihood}
\usage{
ROCnp(t1, d1, t2, d2, b0, t0)
}
\arguments{
    \item{t1}{a vector of length n. Observed times, may be right censored.}
    \item{d1}{a vector of length n, censoring status.
              d=1 means t is uncensored; d=0 means t is right censored. }
    \item{t2}{a vector of length m. Observed times, may be right censored.}
    \item{d2}{a vector of length m, censoring status.}
    \item{b0}{a scalar between 0 and 1. }
    \item{t0}{a scalar, betwenn 0 and 1. }
}
\description{
    Use the empirical likelihood ratio to test the
    hypothesis Ho: (1-b0)th quantile of sample 1 = (1-t0)th quantile 
    of sample 2.
    This is the same as testing Ho: R(t0)= b0, where R(.) is the ROC curve.

The log empirical likelihood been maximized is
\deqn{ \sum_{d1=1} \log \Delta F_1(t1_i) + \sum_{d1=0} \log [1-F_1(t1_i)] 
     + \sum_{d2=1} \log \Delta F_2(t2_j) + \sum_{d2=0} \log [1-F_2(t2_j)] .}
}
\details{
Basically, we first obtain two log likelihood ratios from two samples, 
in testing a common quantile c. And then we minimize the sum of the
two test statistic over c.
 
See the tech report below.
}
\value{
    A list with the following components:
    \item{"-2LLR"}{the -2 loglikelihood ratio; have approximate chisq 
                  distribution under \eqn{H_o}.}
    \item{cstar}{the estimated common quantile.}
}
\references{
    Zhou, M. and Liang, H (2008). 
       Empirical Likelihood for Hybrid Two Sample Problem with
       Censored Data. Tech. Report. 
}
\author{ Mai Zhou. }
\examples{
#### An example of testing the equality of two medians. No censoring.
ROCnp(t1=rexp(100), d1=rep(1,100), t2=rexp(120), 
                                d2=rep(1,120), b0=0.5, t0=0.5)
}
\keyword{nonparametric}
\keyword{htest}
