\name{bjtestII}
\alias{bjtestII}
\title{Alternative test of the Buckley-James estimator by Empirical Likelihood}
\usage{
bjtestII(y, d, x, beta)
}
\arguments{
    \item{y}{a vector of length N, containing the censored responses.}
    \item{d}{a vector of length N. Either 1's or 0's. 
              d=1 means y is uncensored;
	      d=0 means y is right censored. }
    \item{x}{a matrix of size N by q. }
    \item{beta}{a vector of length q. The value of the regression 
           coefficient to be tested in the model 
            \eqn{Y_i = \beta x_i  + \epsilon_i} }
}
\description{
    Use the empirical likelihood ratio (alternative form) and Wilks theorem to test if the
    regression coefficient is equal to beta, based on the estimating equations.

The log empirical likelihood been maximized is
\deqn{ \sum_{j=1}^n \log p_j ; ~  \sum p_j =1 }
where the probability \eqn{p_j} is for the jth martingale differences of the estimating equations.
}
\details{

The above likelihood should be understood as the likelihood of the martingale difference terms. For the definition of the Buckley-James martingale or estimating equation, please see the (2015) book
in the reference list.

The estimation equations used when maximize the 
empirical likelihood is 
\deqn{ 0 = \sum d_i \Delta F(e_i) (x \cdot m[,i])/(n w_i) }
where \eqn{e_i} is the residuals, other details are described in the reference book of 2015 below.

The final test is carried out by \code{el.test}. So the output is similar to the output of \code{el.test}.
}
\value{
    A list with the following components:
    \item{"-2LLR"}{the -2 loglikelihood ratio; have approximate chisq 
                  distribution under \eqn{H_o}.}
    \item{logel2}{the log empirical likelihood, under estimating equation.}
    \item{logel}{the log empirical likelihood of the Kaplan-Meier of e's.}
    \item{prob}{the probabilities that max the empirical likelihood 
               under estimating equation.}
}
\references{
    Zhou, M. (2016) Empirical Likelihood Methods in Survival
        Analysis. CRC Press.
 
    Buckley, J. and James, I. (1979). Linear regression with censored data.
      Biometrika,  \bold{66} 429-36.

    Zhou, M. and Li, G. (2008). Empirical likelihood analysis 
      of the Buckley-James estimator.
      Journal of Multivariate Analysis, \bold{99},  649--664.
    
    Zhu, H. (2007) Smoothed Empirical Likelihood for Quantiles
  and Some Variations/Extension of Empirical Likelihood for Buckley-James Estimator,
     Ph.D. dissertation, University of Kentucky. 
}
\author{ Mai Zhou. }
\examples{
data(myeloma)
bjtestII(y=myeloma[,1], d=myeloma[,2], x=cbind(1, myeloma[,3]), beta=c(37, -3.4))
}
\keyword{nonparametric}
\keyword{htest}
