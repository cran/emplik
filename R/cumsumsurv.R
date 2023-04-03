cumsumsurv <- function(x){
 ## if (sum(is.na(x))>0) stop('NaNs');  3/2015 MZ
 if(any(is.na(x))) {
  stop('NaNs')
  }
  s = x
  return(.C('cumsumsurv',
     x = as.numeric(x),
	 s = as.numeric(s),
	 LLL = length(x),
	 package = 'emplik')$s)
}