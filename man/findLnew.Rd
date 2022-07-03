\name{findLnew}
\alias{findLnew}
\title{Find the Wilks Confidence Interval Lower Bound from the Given (empirical) Likelihood Ratio Function}
\usage{
findLnew(step=0.003, initStep=0, fun, MLE, level=3.84146, tol=.Machine$double.eps^0.5,...)
}
\description{
This function is similar to \code{findUL2} but here the seeking of lower and upper bound are
separate (the other half is \code{findUnew}).

The reason for this is: sometime we need to supply the \code{fun}
with different nuisance parameter(s) values when seeking Lower or Upper bound.
For example \code{fun} returns the -2LLR for a given parameter of interest, 
but there are additional nuisance parameter need to be profiled out, and we 
need to give a range of the nuisance parameter to be max/min over. This range can be 
very different for parameter near Upper bound vs near Lower bound. In the \code{findUL2},
you have to supply a range really wide that (hopefully) works for both Upper and Lower
bound. Here, with separate \code{findLnew} and \code{findUnew} we can tailor the range
for one end of the confidence interval.

Those nuisance parameter(s)
are supplied via the ... input of this function.

Another improvement is that we used the "extendInt" option of the \code{uniroot}.
So now we can and did used a smaller default step input, compare to \code{findUL2}.
 
This program uses uniroot( ) to find (only) the lower (Wilks) confidence
limit based on the -2 log likelihood ratio, which the required 
input \code{fun} is supposed to supply.


Basically, starting from \code{MLE}, we search on lower
direction, by \code{step} away
from \code{MLE}, until we find values that have -2LLR = level.
(the value of -2LLR at MLE is supposed to be zero.)

At curruent implimentation, only handles one dimesional parameter, i.e. only
confidence intervals, not confidence regions.
}
\arguments{
    \item{step}{a positive number. The starting step size of the search. Reasonable value should be about 1/5 of the SD of MLE. }
    \item{initStep}{a nonnegative number. The first step size of the search. Sometimes, you may want to put a larger innitStep to speed the search.}
    \item{fun}{a function that returns a list. One of the item in the list should be "-2LLR", which is the -2 log (empirical) likelihood ratio. 
 The first input of \code{fun} must be the parameter for which we are seeking the confidence interval. (The MLE or NPMLE of this parameter should be supplied as in the input MLE). The rest of the input to \code{fun} are typically the data.
  If the first input of \code{fun} is set to MLE, then the returned -2LLR should be 0.} 
    \item{MLE}{The MLE of the parameter. No need to be exact, as long as it is inside the confidence interval.}
    \item{level}{an optional positive number, controls the confidence level. Default to 3.84146 = chisq(0.95, df=1). 
	Change to 2.70=chisq(0.90, df=1) to get a 90\% confidence interval. }
   \item{tol}{Error bound of the final result.}
   \item{...}{additional arguments, if any, to pass to \code{fun}.}
}
\value{
    A list with the following components:
    \item{Low}{the lower limit of the confidence interval.}
    \item{FstepL}{the final step size when search lower limit. An indication of the precision.}
    \item{Lvalue}{The -2LLR value of the final \code{Low} value. Should be approximately equal to level. If larger than level, 
	than the confidence interval limit \code{Low} is wrong.}
}
\details{
Basically we repeatedly testing the value of the parameter, until we find those
which the -2 log likelihood value is equal to 3.84146 (or other level, if set differently).
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
findLnew(step=0.1, fun=myfun6, MLE=4.0659, x=x, d=d)
}
\keyword{nonparametric}
\keyword{htest}

