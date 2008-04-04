WRegTest <- function(x, y, delta, beta0, psifun=function(t){t}) {
# The test included is the Empirical Likelihood Ratio test
# for the case weighted estimator of the censored AFT model. 
#
# This estimator assumes a random design or a correlation model as 
# defined by Freedman (1981). 
# This coresponds to the bootstrapping the cases in the linear model.

# This estimator in the AFT model is studied by 
# Koul-Susarla-Van Ryzin (1982), Zhou (1992) and Stute (1993) etc.
# For empirical likelihood analysis, please see Zhou, Bathke, Kim (2006).
# The EL is defined as
# The constraint equation is

# Input:
# x is a matrix of N rows (covariates).
# y is the observed (censored) responses --- a vector of length N.
# delta is a vector of length N. delta =1 means (y) is not censored.
#           delta = 0 means y is right censored, i.e. the true
#        response is larger than y.
#
# Output:
# the -2log ELratio, and the P-value 

n <- length(y)
m <- length(beta0)
xx <- as.matrix(x)
xdim <- dim(xx)
if ( xdim[1] != n ) stop("check dim of x")
if ( m != xdim[2] ) stop("check dim of x and beta0")
if ( length(delta) != n ) stop("check length of delta")

###### define the estimating/constraint function ####

myfun <- function(y, xmat, beta) {
      temp1 <- psifun( y -  as.vector( xmat %*% beta ) )
      return( temp1 * xmat )
      }
##### now test if the estimating function is/are zero ####

temp2 <- el.cen.EM2(x=y,d=delta, fun=myfun, mu=rep(0,m), xmat=xx, beta=beta0)
EL <- temp2$"-2LLR"
return(EL)
}
