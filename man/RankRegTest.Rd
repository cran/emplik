\name{RankRegTest}
\alias{RankRegTest}
\title{Test the AFT model Rank Regression estimator by Empirical Likelihood}
\usage{
RankRegTest(y, d, x, beta, type="Gehan")
}
\arguments{
    \item{y}{a vector of length N, containing the censored responses.}
    \item{d}{a vector (length N) of either 1's or 0's. 
              d=1 means y is uncensored;
	      d=0 means y is right censored. }
    \item{x}{a matrix of size N by q. }
    \item{beta}{a vector of length q. the value of the regression 
               coefficient to be tested in the model 
               \eqn{y_i = \beta x_i  + \epsilon_i} }.
    \item{type}{default to Gehan type. 
               The other option is Logrank type.}
}
\description{
    Use the empirical likelihood ratio and Wilks theorem to test if the
    regression coefficient is equal to beta, based on the rank estimator
    for the AFT model.

The log empirical likelihood been maximized is
\deqn{ \sum_{d=1} \log \Delta F(e_i) + \sum_{d=0} \log [1-F(e_i)];}
where \eqn{e_i} are the residuals.
}
\details{

The above likelihood should be understood as the likelihood of the 
error term, so in the regression model the error epsilon should be iid.


The estimation equation used when maximize the 
empirical likelihood is 
\deqn{ 0 = \sum_i \phi (e_i)  d_i \Delta F(e_i) (x_i - \bar x_i  )/(n w_i) }
which was descrided in detail in the reference below.
}
\value{
    A list with the following components:
    \item{"-2LLR"}{the -2 loglikelihood ratio; should have approximate chisq 
                  distribution under \eqn{H_o}.}
    \item{logel2}{the log empirical likelihood, under estimating equation.}
    \item{logel}{the log empirical likelihood of the Kaplan-Meier of e's.}
    \item{prob}{the probabilities that max the empirical likelihood 
               under rank estimating equation.}
}
\references{
Kalbfleisch, J. and Prentice, R. (2002)
{\em The Statistical Analysis of Failure Time Data}. 2nd Ed.
Wiley, New York.  (Chapter 7)

Jin, Z., Lin, D.Y., Wei, L. J. and Ying, Z. (2003).
Rank-based inference for the accelerated failure time model.
{\em Biometrika}, {\bf 90}, 341-353.

Zhou, M. (2005). Empirical likelihood analysis of the rank
        estimator for the censored accelerated failure time model. 
      Biometrika,  \bold{92}, 492-98.
}
\author{ Mai Zhou. }
\examples{
xx <- c(28,-44,29,30,26,27,22,23,33,16,24,29,24,40,21,31,34,-2,25,19)
}
\keyword{nonparametric}
\keyword{htest}
