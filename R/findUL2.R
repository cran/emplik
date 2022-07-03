findUL2 <- function(step = 0.01, initStep = 0, fun, MLE, level = 3.84146, 
                  tol = .Machine$double.eps^0.5, ...) 
{
    value <- 0
    step1 <- step
    Lbeta <- MLE - initStep
	
	    while (value < level) {
            Lbeta <- Lbeta - step1
            value <- fun(Lbeta, ...)$"-2LLR"
        }
        Lbeta0 <- Lbeta
        Lbeta1 <- Lbeta + step1
        tempfun <- function(beta){
                   return( level - fun(beta, ...)$"-2LLR" )
                      }
        temp1 <- uniroot(tempfun, lower=Lbeta0, upper=Lbeta1, tol=tol)
        Lbeta <- temp1$root
        value1 <- level - temp1$f.root
		
    value <- 0
    Ubeta <- MLE + initStep
    
        while (value < level) {
            Ubeta <- Ubeta + step
            value <- fun(Ubeta, ...)$"-2LLR"
        }
        Ubeta0 <- Ubeta
        Ubeta1 <- Ubeta - step    
        temp2 <- uniroot(tempfun, lower=Ubeta1, upper=Ubeta0, tol=tol)
        Ubeta <- temp2$root
        value <- level - temp2$f.root
  
  return(list(Low = Lbeta, Up = Ubeta, FstepL=temp1$estim.prec, FstepU =temp2$estim.prec, 
        Lvalue = value1, Uvalue = value))
}