\name{ROCnp2}
\alias{ROCnp2}
\title{Test the ROC curve by Empirical Likelihood}
\usage{
ROCnp2(t1, d1, t2, d2, b0, t0)
}
\arguments{
    \item{t1}{a vector of length n. Observed times, sample 1. may be right censored.}
    \item{d1}{a vector of length n, censoring status.
              d=1 means t is uncensored; d=0 means t is right censored. }
    \item{t2}{a vector of length m. Observed times, sample 2. may be right censored.}
    \item{d2}{a vector of length m, censoring status.}
    \item{b0}{a scalar, between 0 and 1. }
    \item{t0}{a scalar, betwenn 0 and 1. }
}
\description{
    Use empirical likelihood ratio to test the
    hypothesis Ho: (1-b0)th quantile of sample 1 = (1-t0)th quantile 
    of sample 2.
    This is the same as testing Ho: R(t0)= b0, where R(.) is the ROC curve.

The log empirical likelihood been maximized is
\deqn{ \sum_{d1=1} \log \Delta F_1(t1_i) + \sum_{d1=0} \log [1-F_1(t1_i)] 
     + \sum_{d2=1} \log \Delta F_2(t2_j) + \sum_{d2=0} \log [1-F_2(t2_j)] .}
 
This empirical likelihood ratio has a chi square limit under Ho.
}
\details{
First, we test (1-b0)th quantile of sample 1 = c
and also test (1-t0)th quantile of sample 2 = c. 
This way we obtain two log likelihood ratios. 

Then we minimize the sum of the
two log likelihood ratios over c.
 
This version use an exhaust search for the minimum (over c). 
Since the objective (log lik) are piecewise constants, the
optimum( ) function in R do not work well.
See the tech report below for details on a similar setting.
}
\value{
    A list with the following components:
    \item{"-2LLR"}{the -2 loglikelihood ratio; have approximate chisq 
                  distribution under \eqn{H_o}.}
    \item{cstar}{the estimated common quantile.}
}
\references{
    Su, Haiyan; Zhou, Mai and Liang, Hua (2011). 
       Semi-parametric Hybrid Empirical Likelihood Inference for  Two Sample Comparison with Censored Data.  
\emph{Lifetime Data Analysis}, \bold{17}, 533-551.
}
\author{Mai Zhou}
\examples{
#### An example of testing the equality of two medians. 
#### No censoring. 
# ROCnp2(t1=rexp(100), d1=rep(1,100), t2=rexp(120), 
#                              d2=rep(1,120), b0=0.5, t0=0.5)
###############################################################
#### This example do not work on the Solaris Sparc machine. 
#### But works fine on other platforms.
###########
#### Next, an example of finding 90 percent confidence 
#### interval of R(0.5)
####  Note: We are finding confidence interval for R(0.5). 
#### So we are testing  
####  R(0.5)= 0.35, 0.36, 0.37, 0.38, etc. try to find 
#### values so that testing R(0.5) = L , U has p-value 
#### of 0.10, then [L,  U] is the 90 percent
####  confidence interval for R(0.5).
#set.seed(123)
#t1 <- rexp(200)
#t2 <- rexp(200)
#ROCnp( t1=t1, d1=rep(1, 200), t2=t2, d2=rep(1, 200), 
#                                b0=0.5, t0=0.5)$"-2LLR"
#### since the -2LLR value is less than  
#### 2.705543 = qchisq(0.9, df=1),  so the 
#### confidence interval contains 0.5.
#gridpoints <- 35:65/100
#ELvalues <- gridpoints
#for(i in 1:31) ELvalues[i] <- ROCnp2(t1=t1, d1=rep(1, 200), 
#      t2=t2, d2=rep(1, 200), b0=gridpoints[i], t0=0.5)$"-2LLR"
#myfun1 <- approxfun(x=gridpoints, y=ELvalues)
#uniroot(f=function(x){myfun1(x)-2.705543}, 
#                                interval= c(0.35, 0.5) )
#uniroot(f= function(x){myfun1(x)-2.705543}, 
#                                interval= c(0.5, 0.65) )
#### So, taking the two roots, we see the 90 percent 
#### confidence interval for R(0.5) in this 
#### case is [0.4457862,  0.5907723].
###############################################
}
\keyword{nonparametric}
\keyword{htest}
