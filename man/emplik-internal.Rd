\name{emplik-internal}
\alias{Wdataclean2}
\alias{DnR}
\alias{logelr}
\alias{llog}
\alias{llogp}
\alias{llogpp}
\alias{solve3.QP}
\alias{WKM}
\alias{el.test.wt}
\title{Internal emplik functions}
\usage{
logelr(x,mu,lam)
llog(z, eps)
llogp(z, eps)
llogpp(z, eps)
Wdataclean2(z,d,wt=rep(1,length(d)))
DnR(x,d,w)
solve3.QP(D, d, A, b, meq, factorized=FALSE)
WKM(x,d,w=rep(1,length(d)))
el.test.wt(x,wt,mu)
}
\description{ Internal emplik functions }
\details{
These are not to be called by the user. 

\code{Wdataclean2} and
\code{DnR} are used by the functions \code{emplikH1.test}, 
\code{emplikH2.test} and \code{emplikdisc.test}.  

\code{logelr}, \code{llog}, \code{llogp} and \code{llogpp} 
are used by function \code{el.test}.

\code{Wdataclean2}, \code{WKM} and \code{solve3.QP} are
used by function \code{el.cen.test}. 

\code{el.test.wt}, \code{Wdataclean2} are used by \code{el.cen.EM}.

\code{el.test.wt} is basically a weighted version of \code{el.test}
but only for 1 dim mean. 
}
\keyword{internal}
