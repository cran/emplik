\name{el.test.wt}
\alias{el.test.wt}
\title{Weighted Empirical Likelihood ratio for mean, uncensored data}
\usage{
el.test.wt(x, wt, mu)
}
\description{
This program is similar to \code{el.test( )} except
it takes weights, and is for one dimensional mu.

The mean constraint considered is:
\deqn{ \sum_{i=1}^n  p_i x_i  = \mu . }
where \eqn{p_i = \Delta F(x_i)} is a probability.
Plus the probability constraint: \eqn{ \sum p_i =1}.

The weighted log empirical likelihood been maximized is
\deqn{ \sum_{i=1}^n w_i \log p_i. } 
}
\arguments{
    \item{x}{a vector containing the observations.}
    \item{wt}{a vector containing the weights.}
    \item{mu}{a real number used in the constraint, weighted
                mean value of \eqn{f(X)}.}
}
\value{
    A list with the following components:
    \item{x}{the observations.}
    \item{wt}{the vector of weights.}
    \item{prob}{The probabilities that maximized the weighted 
                  empirical likelihood under mean constraint.}
}
\details{
This function used to be an internal function. 
It becomes external because others may find it useful
elsewhere.

The constant \code{mu} must be inside 
\eqn{( \min x_i , \max x_i ) }
for the computation to continue. 

}
\author{ Mai Zhou }
\references{
    Zhou, M. (2002). 
        Computing censored empirical likelihood ratio 
        by EM algorithm. 
    \emph{Tech Report, Univ. of Kentucky, Dept of Statistics}
}
\examples{
## example with tied observations
x <- c(1, 1.5, 2, 3, 4, 5, 6, 5, 4, 1, 2, 4.5)
d <- c(1,   1, 0, 1, 0, 1, 1, 1, 1, 0, 0,   1)
el.cen.EM(x,d,mu=3.5)
## we should get "-2LLR" = 1.2466....
myfun5 <- function(x, theta, eps) {
u <- (x-theta)*sqrt(5)/eps 
INDE <- (u < sqrt(5)) & (u > -sqrt(5)) 
u[u >= sqrt(5)] <- 0 
u[u <= -sqrt(5)] <- 1 
y <- 0.5 - (u - (u)^3/15)*3/(4*sqrt(5)) 
u[ INDE ] <- y[ INDE ] 
return(u)
}
el.cen.EM(x, d, fun=myfun5, mu=0.5, theta=3.5, eps=0.1)
}
\keyword{nonparametric}
\keyword{htest}
