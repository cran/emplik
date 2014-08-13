cumsumsurv<-function(x){
 	if (sum(is.na(x))>0)stop('NaNs');
 	s=x;
 	.C('cumsumsurv',x=as.numeric(x),s=as.numeric(s),LLL=length(x))$s
 	}



