#####################################################################
## Improved by Mai Zhou 10/10/2021 by using "extendInt" of uniroot( ).
## The result is that it is not sensitive to step value (in terms of computing time). 
## We can use a smaller step= value without slowing down much.
## Also, Upper and Lower bound seeking are now separate. Since some parameter(s) 
## in ... for fun can take different values when seeking Upper/Lower bound. An 
## example of this is the pAUC estimation.
#####################################################################
findUnew <- function(step=0.003, initStep=0, fun, MLE, level=3.84146, tol=.Machine$double.eps^0.5, ...)
{
Ubeta <- MLE + initStep
Ubeta0 <- Ubeta + 3*step
Ubeta1 <- Ubeta

tempfun <- function(beta){return(level-fun(beta,...)$"-2LLR")}
temp2 <- uniroot(tempfun,lower=Ubeta1,upper=Ubeta0, extendInt="downX", tol=tol)
Ubeta <- temp2$root
value <- level-temp2$f.root
return(list(Up=Ubeta, FstepU=temp2$estim.prec, Uvalue=value))
}
