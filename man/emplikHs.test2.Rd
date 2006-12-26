\name{emplikHs.test2}
\alias{emplikHs.test2}
\title{Two sample empirical likelihood ratio test for hazards
with right censored, left truncated data. Many constraints.}
\usage{
emplikHs.test2(x1, d1, y1= -Inf, x2, d2, y2 = -Inf,
          theta, fun1, fun2, maxit=25,tola = 1e-7,itertrace =FALSE)
}
\description{
Use empirical likelihood ratio and Wilks theorem to test 
the null hypothesis that
\deqn{ 
\int{f_1(t) dH_1(t))} -
\int{f_2(t) dH_2(t))} = \theta 
}
where \eqn{H_*(t)} is the (unknown) cumulative
hazard functions; \eqn{f_*(t)} can be any predictable  
functions of \eqn{t}. 
\eqn{\theta} is a vector of parameters (dim=q). 
The given value of \eqn{\theta}
in these computation are the value to be tested.
The data can be right censored and left truncated.

When the given constants \eqn{\theta} is too far
away from the NPMLE, there will be no hazard function satisfy this 
constraint and the -2 Log empirical likelihood ratio
will be infinite. In this case the computation will stop.
}
\arguments{
    \item{x1}{a vector, the observed survival times, sample 1.}
    \item{d1}{a vector, the censoring indicators, 1-uncensor; 0-censor.}
    \item{y1}{optional vector, the left truncation times.}
    \item{x2}{a vector, the observed survival times, sample 2.}
    \item{d2}{a vector, the censoring indicators, 1-uncensor; 0-censor.}
    \item{y2}{optional vector, the left truncation times.}
    \item{fun1}{a predictable function used to calculate
         the weighted discrete hazard in \eqn{H_0}. 
         \code{fun1(x)} must be able to take a vector input (length n)
         \code{x}, and output a matrix of n x q.}
    \item{fun2}{ Ditto.}
    \item{tola}{an optional positive real number, the tolerance of
       iteration error in solve the non-linear equation needed in constrained 
        maximization.}
    \item{theta}{a given vector of length q. for Ho constraint. }
    \item{maxit}{integer, maximum number of iteration. }
    \item{itertrace}{ Logocal, lower bound for lambda }
}
\value{
    A list with the following components:
    \item{"-2LLR"}{The -2Log Likelihood ratio.}
    \item{lambda}{the final value of the Lagrange multiplier.}
    \item{"-2LLR(sample1)"}{The -2Log empirical likelihood ratio for sample 
                            one only. Useful in one sample problems.}
    \item{niters}{number of iterations used}
}
\details{
The log likelihood been maximized is the Poisson likelihood:
\deqn{ 
\sum D_{1i} \log w_i  - \sum R_{1i} w_i + 
\sum D_{2j} \log v_j  - \sum R_{2j} v_j 
}
where \eqn{w_i = \Delta H_1(t_i)} is the jump 
of the cumulative hazard function at \eqn{t_i}
(for first sample), 
\eqn{D_{1i}} is the number of failures
observed at \eqn{t_i}, \eqn{R_{1i}} is 
the number of subjects at risk at
time \eqn{t_i}. Dido for sample two.

For (proper)
discrete distributions, the jump size of the cumulative hazard at
the last jump is always 1. So, in the likelihood ratio, it cancels.
But the last jump of size 1 still matter when computing the constraint.

The constants \code{theta} must be inside the so called
feasible region for the computation to continue. This is similar to the
requirement that in testing the value of the mean, the value must be
inside the convex hull of the observations.
It is always true that the NPMLE values are feasible. So when the
computation stops, try move the \code{theta} closer
to the NPMLE, which we print out first thing in this function,
even when other later computations do not go. 
 When the computation stops, the -2LLR should have value
infinite.

You can also use this function for one sample problems. 
You need to artificially supply data for sample two of minimal size
(like size 2q+2), and specify a fun2() that ALWAYS return 0
(zero vector, or zero matrix).
Then, look for -2LLR(sample1) in the output.

}
\author{ Mai Zhou }
\references{
    Zhou and Fang (2001). 
    ``Empirical likelihood ratio for 2 sample problems for censored data''. 
    \emph{Tech Report, Univ. of Kentucky, Dept of Statistics}
}
\examples{
if(require("boot", quietly = TRUE)) {
####library(boot)
data(channing)
ymale <- channing[1:97,2]
dmale <- channing[1:97,5]
xmale <- channing[1:97,3]
yfemale <- channing[98:462,2]
dfemale <- channing[98:462,5]
xfemale <- channing[98:462,3]
fun1 <- function(x) { as.numeric(x <= 960) }
########################################################
fun2 <- function(x){ cbind(as.numeric(x <= 960), as.numeric(x <= 860))}
############ fun2 has matrix output ###############
emplikHs.test2(x1=xfemale, d1=dfemale, y1=yfemale, x2=xmale,d2=dmale, y2=ymale, theta=c(0,0), fun1=fun2, fun2=fun2)
}
#############################################
###################### Second example:
if(require("KMsurv", quietly = TRUE)) {
####library(KMsurv)
data(kidney)
### these functions counts the risk set size, so delta=1 always ###
temp1 <- Wdataclean3(z=kidney$time[kidney[,3]==1], d=rep(1,43) )
temp2 <- DnR(x=temp1$value, d=temp1$dd, w=temp1$weight)
TIME <- temp2$times
RISK <- temp2$n.risk
fR1 <- approxfun(x=TIME, y=RISK, method="constant", yright=0, rule=2, f=1)
temp1 <- Wdataclean3(z=kidney$time[kidney[,3]==2], d=rep(1,76) )
temp2 <- DnR(x=temp1$value, d=temp1$dd, w=temp1$weight)
TIME <- temp2$times
RISK <- temp2$n.risk
fR2 <- approxfun(x=TIME, y=RISK, method="constant", yright=0, rule=2, f=1)

### the weight function for two sample Gehan-Wilcoxon type test ###
fun <- function(t){ fR1(t)*fR2(t)/((76*43)*sqrt(119/(76*43)) )}
### Here comes the test: ###
emplikHs.test2(x1=kidney[kidney[,3]==1,1],d1=kidney[kidney[,3]==1,2],
   x2=kidney[kidney[,3]==2,1],d2=kidney[kidney[,3]==2,2],
   theta=0, fun1= fun, fun2=fun)
### The results should include this ###
#$"-2LLR"
#[1] 0.002473070
#
#$lambda
#[1] -0.1713749
#######################################
######### the weight function for log-rank test #####
funlogrank <- function(t){sqrt(119/(76*43))*fR1(t)*fR2(t)/(fR1(t)+fR2(t))}
##### Now the log-rank test ###
emplikHs.test2(x1=kidney[kidney[,3]==1,1],d1=kidney[kidney[,3]==1,2],
  x2=kidney[kidney[,3]==2,1],d2=kidney[kidney[,3]==2,2],
  theta=0, fun1=funlogrank, fun2=funlogrank)
##### The result of log rank test should include this ###
#
#$"-2LLR"
#[1] 2.655808
#
#$lambda
#[1] 3.568833
#######################################################
###### the weight function for both type test ####
funBOTH <- function(t) {
       cbind(sqrt(119/(76*43))*fR1(t)*fR2(t)/(fR1(t)+fR2(t)), 
                 fR1(t)*fR2(t)/((76*43)*sqrt(119/(76*43)))) }
#### The test that combine both tests ###
emplikHs.test2(x1=kidney[kidney$type==1,1],d1=kidney[kidney$type==1,2],
    x2=kidney[kidney$type==2,1],d2=kidney[kidney$type==2,2],
    theta=c(0,0), fun1=funBOTH, fun2=funBOTH)
#### the result should include this ###
#
#$"-2LLR"
#[1] 13.25476
#
#$lambda
#[1]  14.80228 -21.86733
##########################################
}
}
\keyword{nonparametric}
\keyword{survival}
\keyword{htest}
