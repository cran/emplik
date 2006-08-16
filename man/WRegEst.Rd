\name{WRegEst}
\alias{WRegEst}
\title{Compute the case weighted regression estimator for AFT model}
\usage{
WRegEst(x, y, delta, LS=TRUE, tau=0.5)
}
\arguments{
    \item{x}{a matrix of size N by q. }
    \item{y}{a vector of length N, containing the censored responses.}
    \item{delta}{a vector (length N) of either 1's or 0's. 
              d=1 means y is uncensored;
	      d=0 means y is right censored. }
    \item{LS}{a logical value. If TRUE then the function will 
              return the least squares estimator. If FALSE then the
              function will return the regression quantile estimator,
              with quantile specified by tau. }.
    \item{tau}{a scalar, between 0 and 1. The quantile to be used in
                 quantile regression. }
}
\description{
For the AFT model, this function computes the case weighted estimator of
beta. Either the least squares estimator or the regression quantile estimator.
}
\details{

The estimator is the minimizer of
\deqn{ \sum_{i=1}^n w_i \rho (Y_i - X_i b) }

Assuming a correlation model
\eqn{ Y_i = X_i \beta + \sigma(X_i) \epsilon_i }.

}
\value{
the estimator \eqn{ \hat \beta}.
}
\references{
    Zhou, M.; Bathke, A. and Kim, M. (2006). 
       Empirical likelihood analysis of the
       Heteroscastic AFT model. Tech. Report. 
}
\author{ Mai Zhou. }
\examples{
#library(rankreg)
#data(smallcell)
#WRegEst(x=cbind(1,smallcell[,1],smallcell[,2]), 
#        y=smallcell[,3], delta=smallcell[,4]) 
####################################################
#### you should get         x1         x2         x3
####                 -59.22126 -488.41306   16.03259
####################################################
xx <- c(28,-44,29,30,26,27,22,23,33,16,24,29,24,40,21,31,34,-2,25,19)
}
\keyword{nonparametric}
