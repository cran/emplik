\name{emplik-internal}
\alias{Wdataclean2}
\alias{Wdataclean3}
\alias{Wdataclean5}
\alias{DnR}
\alias{logelr}
\alias{logwelr}
\alias{gradf}
\alias{llog}
\alias{llogp}
\alias{llogpp}
\alias{solve3.QP}
\alias{WKM}
\alias{el.test.wt}
\alias{el.test.wt2}
\alias{el.test.wt3}
\alias{LTRC}
\title{Internal emplik functions}
\usage{
logelr(x,mu,lam)
logwelr(x,mu,wt,lam)
gradf(z,wt,lam)
llog(z, eps)
llogp(z, eps)
llogpp(z, eps)
Wdataclean2(z,d,wt=rep(1,length(d)))
Wdataclean3(z,d,zc=1:length(z),wt=rep(1,length(d)))
Wdataclean5(z,d,zc=1:length(z),wt=rep(1,length(d)),xmat)
DnR(x,d,w,y=rep(-Inf,length(x)))
solve3.QP(D, d, A, b, meq, factorized=FALSE)
WKM(x,d,zc=rep(1,length(x)),w=rep(1,length(d)))
LTRC(x,d,w=rep(1, length(d)),y=rep(-Inf, length(x)))
el.test.wt(x,wt,mu)
el.test.wt2(x,wt,mu,maxit,gradtol,Hessian,svdtol,itertrace)
el.test.wt3(x,wt,mu,maxit,gradtol,Hessian,svdtol,itertrace)
}
\description{ Internal emplik functions }
\details{
These are not to be called by the user. 

\code{Wdataclean2} and
\code{DnR} are used by the functions \code{emplikH1.test}, 
\code{emplikH2.test} and \code{emplikdisc.test}.  It is also used by 
\code{LTRC}.  

\code{logelr}, \code{llog}, \code{llogp} and \code{llogpp} 
are used by function \code{el.test}.

\code{Wdataclean2}, \code{WKM} and \code{solve3.QP} are
used by function \code{el.cen.test}. 

\code{el.test.wt}, \code{Wdataclean2} are used by \code{el.cen.EM}.

\code{el.test.wt} is basically a weighted version of \code{el.test}
but only for 1 dim mean. 

\code{LTRC} is for Left Truncated and Right Censored data.

\code{el.test.wt2}, \code{el.test.wt3} are similar to \code{el.test.wt} 
but can take vector mean as constraint. \code{llog}, \code{llogp} and
\code{llogpp} are used by both functions. In addition \code{logwelr}
is used by \code{el.test.wt3}.

}
\keyword{internal}
