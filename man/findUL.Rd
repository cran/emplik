\name{findUL}
\alias{findUL}
\title{Find the Wilks Confidence Interval from the Given (empirical) Likelihood Ratio Function}
\usage{
findUL (step = 0.01, initStep =0, fun, MLE, level = 3.84146, ...)
}
\description{
This program uses uniroot( ) to find the upper and lower (Wilks) confidence
limits based on the -2 log likelihood ratio, which the required input \code{fun} is supposed to supply.

Basically, starting from \code{MLE}, we search on both 
directions, by \code{step} away
from \code{MLE}, until we find values that have -2LLR = level.
(the value of -2LLR at MLE is supposed to be zero.)

At curruent implimentation, only handles one dimesional parameter, i.e. only
confidence intervals, not confidence regions.

For examples of using this function to find confidence interval,
see the pdf vignettes file.
}
\arguments{
    \item{step}{a positive number. The starting step size of the search. Reasonable value should be about 1/5 of the SD of MLE. }
    \item{initStep}{a nonnegative number. The first step size of the search. Sometimes, you may want to put a larger innitStep to speed the search.}
    \item{fun}{a function that returns a list. One of the item in the list should be "-2LLR", which is the -2 log (empirical) likelihood ratio. 
 The first input of \code{fun} must be the parameter for which we are seeking the confidence interval. (The MLE or NPMLE of this parameter should be supplied as in the input MLE). The rest of the input to \code{fun} are typically the data.
  If the first input of \code{fun} is set to MLE, then the returned -2LLR should be 0.} 
    \item{MLE}{The MLE of the parameter. No need to be exact, as long as it is inside the confidence interval.}
    \item{level}{an optional positive number, controls the confidence level. Default to 3.84 = chisq(0.95, df=1). Change to 2.70=chisq(0.90, df=1) to get a 90\% confidence interval. }
   \item{...}{additional arguments, if any, to pass to \code{fun}.}
}
\value{
    A list with the following components:
    \item{Low}{the lower limit of the confidence interval.}
    \item{Up}{the upper limit of the confidence interval.}
    \item{FstepL}{the final step size when search lower limit. An indication of the precision.}
    \item{FstepU}{Ditto. An indication of the precision of the upper limit.}
    \item{Lvalue}{The -2LLR value of the final \code{Low} value. Should be approximately equal to level. If larger than level, than the confidence interval limit \code{Low} is wrong.}
    \item{Uvalue}{Ditto. Should be approximately equal to \code{level}.}
}
\details{

Basically we repeatedly testing the value of the parameter, until we find those
which the -2 log likelihood value is equal to 3.84 (or other level, if set differently).

If there is no value exactly equal to 3.84, it is better to use findULold( ), where
we stop at the value which result
a -2 log likelihood just below 3.84. (as in the discrete case, like quantiles.)

}
\author{ Mai Zhou }
\references{
    Zhou, M. (2016). Empirical Likelihood Method in Survival Analysis. CRC Press.
}
\examples{
## example with tied observations. Kaplan-Meier mean=4.0659.
## For more examples see vignettes.
x <- c(1, 1.5, 2, 3, 4, 5, 6, 5, 4, 1, 2, 4.5)
d <- c(1,   1, 0, 1, 0, 1, 1, 1, 1, 0, 0,  1)
myfun6 <- function(theta, x, d) {
el.cen.EM2(x, d, fun=function(t){t}, mu=theta)
}
findUL(step=0.2, fun=myfun6, MLE=4.0659, x=x, d=d)
}
\keyword{nonparametric}
\keyword{htest}

