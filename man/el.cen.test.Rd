\name{el.cen.test}
\alias{el.cen.test}
\title{Empirical likelihood ratio for mean
with right censored data, by QP.}
\usage{
el.cen.test(x,d,fun=function(x){x},mu,error=1e-8,maxit=20)
}
\description{
This program computes the maximized (wrt \eqn{p_i}) empirical
log likelihood function for right censored data with 
the MEAN constraint:
\deqn{ \sum_i [ d_i p_i g(x_i) ] = \int g(t) dF(t) = \mu  } 
where \eqn{p_i = \Delta F(x_i)} is a probability,
\eqn{d_i} is the censoring indicator. 
The \eqn{d} for the largest observation is always taken to be 1.
It then computes the -2 log 
empirical likelihood ratio which should be approximately chi-square
distributed if the constraint is true.
Here \eqn{F(t)} is the (unknown) CDF; 
\eqn{g(t)} can be any given left
continuous function in \eqn{t}.
\eqn{\mu} is a given constant. 
The data must contain some right censored observations.
If there is no censoring or the only censoring is the largest 
observation, the code will stop and we should use 
\code{el.test} that is for uncensored data.  

The log likelihood been maximized is
\deqn{ \sum_{d_i=1} \log \Delta F(x_i) + \sum_{d_i=0} \log [ 1-F(x_i) ].}
}
\arguments{
    \item{x}{a vector containing the observed survival times.}
    \item{d}{a vector containing the censoring indicators, 
           1-uncensor; 0-censor.}
    \item{fun}{a left continuous (weight) function used to calculate
         the mean as in \eqn{H_0}. 
         \code{fun(t)} must be able to take a vector input \code{t}. 
      Default to the identity function \eqn{f(t)=t}.}
    \item{mu}{a real number used in the constraint, sum to this value.}
    \item{error}{an optional positive real number specifying the tolerance of
       iteration error in the QP. This is the bound of the
       \eqn{L_1} norm of the difference of two successive weights.}
    \item{maxit}{an optional integer, used to control maximum number of
             iterations. }
}
\value{
    A list with the following components:
    \item{"-2LLR"}{The -2Log Likelihood ratio.}
    \item{xtimes}{the location of the hazard jumps.}
    \item{weights}{the jump size of CDF function at those locations.}
    \item{Pval}{P-value}
    \item{error}{the \eqn{L_1} norm between the last two \code{wts}.}
    \item{iteration}{number of iterations carried out}
}
\details{
When the given constants \eqn{\mu} is too far
away from the NPMLE, there will be no distribution
satisfy the constraint.
In this case the computation will stop.
The -2 Log empirical likelihood ratio
should be infinite. 

The constant \code{mu} must be inside 
\eqn{( \min f(x_i) , \max f(x_i) ) }
for the computation to continue. 
It is always true that the NPMLE values are feasible. So when the
computation stops, try move the \code{mu} closer
to the NPMLE, or use a different \code{fun}. 

This function uses sequential Quadratic Programming to find the
maximum. Unlike other functions in this package,
it can be slow for larger sample sizes.
It took about one minute for a sample of size 2000 with 20\% censoring
on a 1GHz, 256MB PC, about 19 seconds on a 3 GHz 512MB PC.
}
\author{ Mai Zhou, Kun Chen }
\references{
    Pan, X. and Zhou, M. (1999). Empirical likelihood ratio, one parameter
                         sub-family of distributions and censored data.
    \emph{J. Statist. Plann. Inference}. \bold{75}, 379-392.

    Chen, K. and Zhou, M. (2000). 
        Computing censored empirical likelihood ratio 
        using Quadratic Programming. 
    \emph{Tech Report, Univ. of Kentucky, Dept of Statistics}
}
\examples{
el.cen.test(rexp(100), c(rep(0,25),rep(1,75)), mu=1.5)
## second example with tied observations
x <- c(1, 1.5, 2, 3, 4, 5, 6, 5, 4, 1, 2, 4.5)
d <- c(1,   1, 0, 1, 0, 1, 1, 1, 1, 0, 0,   1)
el.cen.test(x,d,mu=3.5)
# we should get  "-2LLR" = 1.246634  etc. 
}
\keyword{nonparametric}
\keyword{survival}
\keyword{htest}


{  # This program compute the maximized (wrt w_i) loglikelihood 
   # function for right censored data with the MEAN constraint: 
   #                   sum[d_i * w_i * fun(x_i)]= mu        (star)
   # where w_i is the probability.
   # Inputs
   # x  : life time data vector, size n
   # d  : censoring indicator. length n vector. d=1 uncensor; d=0 censor.
   # fun : function to compute the mean, default is f(t)=t. 
   # mu : constrained mean value, a scalar. 
   # error : iteration stopping control standard. default=1e-8.
   # maxit : maximum number of iteration
   #
   # Output
   # -2LLR   : -2log maximized loglikelihood ratio (chisq under Ho)
   # weight  : weight under mean constraint
   # m       : the number of iterations
   #************************************************************
   # depends on Wdataclean2(), WKM() and solve3.QP() 
   ####################################################################
} 
