\name{el.test.wt2}
\alias{el.test.wt2}
\title{Weighted Empirical Likelihood ratio for mean(s), uncensored data}
\usage{
el.test.wt2(x, wt, mu, maxit = 25, gradtol = 1e-07, Hessian = FALSE, 
    svdtol = 1e-09, itertrace = FALSE)
}
\description{
This program is similar to el.test except
it takes weights.

The mean constraint are:
\deqn{ \sum_{i=1}^n  p_i x_i  = \mu . }
Where \eqn{p_i = \Delta F(x_i)} is a probability.
Plus the probability constraint: \eqn{ \sum p_i =1}.

The weighted log empirical likelihood been maximized is
\deqn{ \sum_{i=1}^n w_i \log p_i. } 
}
\arguments{
    \item{x}{a matrix (of size nxp) or vector containing the observations.}
    \item{wt}{a vector of length n, containing the weights.}
    \item{mu}{a vector of length p, used in the constraint. weighted
                mean value of \eqn{f(X)}.}
    \item{maxit}{the maximum number of iteration.}
    \item{gradtol}{the tol for a solution}
    \item{Hessian}{if the Hessian needs to be computed?}
    \item{svdtol}{ tol in perform SVD of the matrix. }
    \item{itertrace}{TRUE/FALSE, if the intermediate steps needs to be printed.}
}
\value{
    A list with the following components:
    \item{lambda}{the Lagrange multiplier.}
    \item{wt}{the vector of weights.}
    \item{prob}{The probabilities that maximized the weighted 
                  empirical likelihood under mean constraint.}
}
\details{
This function used to be an internal function. 
It becomes external because others may find it useful.

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
