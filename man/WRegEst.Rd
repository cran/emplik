\name{WRegEst}
\alias{WRegEst}
\title{Compute the casewise weighted regression estimator for AFT model}
\usage{
WRegEst(x, y, delta, LS=TRUE, tau=0.5)
}
\arguments{
    \item{x}{a matrix of size N by q. }
    \item{y}{a vector of length N, containing the censored responses. Usually 
             the log of the original observed failure times. }
    \item{delta}{a vector (length N) of either 1's or 0's. 
              d=1 means y is uncensored;
	      d=0 means y is right censored. }
    \item{LS}{a logical value. If TRUE then the function will 
              return the least squares estimator. If FALSE then the
              function will return the quantile regression estimator,
              with the quantile level specified by \code{tau}. }
    \item{tau}{a scalar, between 0 and 1. The quantile to be used in
                 quantile regression. If tau=0.5 then it is the median regression. If LS=TRUE, then it is ignored.}
}
\description{
For the AFT model, this function computes the case weighted estimator of
beta. Either the least squares estimator or the regression quantile estimator.
}
\details{
Due to the readily available minimizer, we only provide least squares
and quantile regression here. However, in the companion testing function
\code{WRegTest } the user can supply a self defined \code{psifun} function,
corresponding to the general M-estimation in the regression modeling.
(since there is no minimization needed).

The estimator is the minimizer of
\deqn{ \sum_{i=1}^n w_i \rho (Y_i - X_i b) }
Assuming a correlation model
\deqn{ Y_i = X_i \beta + \sigma(X_i) \epsilon_i }
where \eqn{\rho( )} is either the square or the absolute value function.
}
\value{
The estimator \eqn{ \hat \beta}.
}
\references{
    Zhou, M.; Bathke, A. and Kim, M. (2012). 
       Empirical likelihood analysis of the
       Heteroscastic Accelerated Failure Time model. 
\emph{Statistica Sinica}, \bold{22}, 295-316.
}
\author{ Mai Zhou. }
\examples{
data(smallcell)
WRegEst(x=cbind(1,smallcell[,1],smallcell[,2]), 
        y=smallcell[,3], delta=smallcell[,4]) 
####################################################
#### you should get         x1         x2         x3
####                 -59.22126 -488.41306   16.03259
####################################################
WRegEst(x=cbind(1,smallcell[,1],smallcell[,2]), 
        y=log10(smallcell[,3]), delta=smallcell[,4], LS=FALSE) 
########################################################
#### you should get      
####     [1]     2.603342985  -0.263000044  0.003836832
########################################################
}
\keyword{nonparametric}
