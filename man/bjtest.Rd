\name{bjtest}
\alias{bjtest}
\title{Test the Buckley-James estimator by Empirical Likelihood}
\usage{
bjtest(y, d, x, beta)
}
\arguments{
    \item{y}{a vector of length N, containing the censored responses.}
    \item{d}{a vector (length N) of either 1's or 0's. 
              d=1 means y is uncensored;
	      d=0 means y is right censored. }
    \item{x}{a matrix of size N by q. }
    \item{beta}{a vector of length q. the value of the regression coefficient to 
          be tested in the model \eqn{y_i = \beta x_i  + \epsilon_i} }.
}
\description{
    Use the empirical likelihood ratio and Wilks theorem to test if the
    regression coefficient is equal to beta.

The log empirical likelihood been maximized is
\deqn{ \sum_{d=1} \log \Delta F(e_i) + \sum_{d=0} \log [1-F(e_i)];}
where \eqn{e_i} are the residuals.
}
\details{

The above likelihood should be understood as the likelihood of the 
error term, so in the regression model the error epsilon should be iid.

This version can handle the model where beta is a vector (of length q).

The estimation equation used when maximize the 
empirical likelihood is 
\deqn{ \sum d_i \Delta F(e_i) (x \cdot m[,i])/(n w_i) }
which was descrided in detail in the reference below.
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
    Buckley, J. and James, I. (1979). Linear regression with censored data.
      Biometrika,  \bold{66} 429-36.

    Zhou, M. and Li, G. (2004). Empirical likelihood analysis of the
       Buckley-James estimator. Tech. Report. 
}
\author{ Mai Zhou. }
\examples{
xx <- c(28,-44,29,30,26,27,22,23,33,16,24,29,24,40,21,31,34,-2,25,19)
}
\keyword{nonparametric}
\keyword{htest}
