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

The estimator of beta can be obtained by function 
\code{rankaft( )}
in the package \code{rankreg}. But here you may test other values of 
beta. If you test the beta value that is obtained from the \code{rankaft( )},
then the -2LLR should be 0 and the p-value should be 1.

The above likelihood should be understood as the likelihood of the 
error term, so in the regression model the error \eqn{e_i} should be iid.

The estimation equation used when maximize the 
empirical likelihood is 
\deqn{ 0 = \sum_i \phi (e_i)  d_i \Delta F(e_i) (x_i - \bar x_i  )/(n w_i) }
which was described in detail in the references below.

See also the function \code{RankRegTestH}, which is based on the hazard likelihood.
}
\value{
    A list with the following components:
    \item{"-2LLR"}{the -2 loglikelihood ratio; should have approximate chisq 
                  distribution under \eqn{H_o}.}
    \item{logel2}{the log empirical likelihood, under estimating equation.}
    \item{logel}{the log empirical likelihood of the Kaplan-Meier of e's.}
    \item{prob}{the probabilities that max the empirical likelihood 
               under rank estimating equation constraint.}
}
\references{
Kalbfleisch, J. and Prentice, R. (2002)
The Statistical Analysis of Failure Time Data. 2nd Ed.
Wiley, New York.  (Chapter 7)

Jin, Z., Lin, D.Y., Wei, L. J. and Ying, Z. (2003).
Rank-based inference for the accelerated failure time model.
Biometrika, \bold{90}, 341-53.

Zhou, M. (2005). Empirical likelihood analysis of the rank
        estimator for the censored accelerated failure time model. 
      Biometrika, \bold{92}, 492-98.
}
\author{ Mai Zhou. }
\examples{
data(myeloma)
RankRegTest(y=myeloma[,1], d=myeloma[,2], x=myeloma[,3], beta= -15.50147)
# you should get "-2LLR" = 9.050426e-05 (practically zero)
# The beta value, -15.50147, was obtained by rankaft() from the rankreg package.
}
\keyword{nonparametric}
\keyword{htest}
