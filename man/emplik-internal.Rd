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
\alias{WCY}
\alias{el.test.wt3}
\alias{LTRC}
\alias{iter}
\alias{MMRtime}
\alias{redistF}
\alias{gradf2}
\alias{gradf3}
\title{Internal emplik functions}
\usage{
logelr(x,mu,lam)
logwelr(x,mu,wt,lam)
gradf(z,wt,lam)
llog(z, eps)
llogp(z, eps)
llogpp(z, eps)
Wdataclean2(z,d,wt=rep(1,length(z)))
Wdataclean3(z,d,zc=rep(1, length(z)),wt=rep(1,length(z)))
Wdataclean5(z,d,zc=rep(1, length(z)),wt=rep(1,length(z)),xmat)
DnR(x,d,w,y=rep(-Inf,length(x)))
solve3.QP(D, d, A, b, meq, factorized=FALSE)
WKM(x,d,zc=rep(1,length(d)),w=rep(1,length(d)))
WCY(x,d,zc=rep(1,length(d)),wt=rep(1,length(d)),maxit=25,error=1e-09)
LTRC(x,d,w=rep(1, length(d)),y=rep(-Inf, length(x)))
el.test.wt3(x,wt,mu,maxit,gradtol,Hessian,svdtol,itertrace)
iter(x, y, delta, beta)
MMRtime(x, d, age)
redistF(y, d, Fdist)
gradf2(lam, funt1, evt1, rsk1, funt2, evt2, rsk2, K, n)
gradf3(lam, funt1, evt1, rsk1, funt2, evt2, rsk2, K, n)
}
\description{ Internal emplik functions }
\details{
These are not intended to be called by the user. 

\code{Wdataclean2} and
\code{DnR} are used by the functions \code{emplikH1.test}, 
\code{emplikH2.test} and \code{emplikdisc.test}.  It is also used by 
\code{LTRC}.  

\code{DnR} take the survival data as input, compute the number at risk 
and number of observed events, and output only those time points where
there is some events. Notice the input must be ordered! We should use
Wdataclean2/3/5( ) first and then use the outpt there as the input of DnR.

\code{Wdataclean2} will sort the data and collaps those that are true ties,
and the number of tied value is in the weights.

\code{logelr}, \code{llog}, \code{llogp} and \code{llogpp} 
are used by function \code{el.test}. They are from Owen.

\code{Wdataclean2}, \code{WKM} and \code{solve3.QP} are
used by function \code{el.cen.test}. 
WKM() is the weighted Kaplan-Meier.

\code{WCY} is the weighted Chang and Yang self-consistent estimator
for doubly censored data.

\code{Wdataclean2} is used by \code{el.cen.EM}.

\code{LTRC} is for Left Truncated and Right Censored data.

\code{gradf3} is used by the function emplikHs.test2. (and stored there)

\code{el.test.wt3} are similar to \code{el.test.wt} and \code{el.test.wt2}, 
but can take vector mean as constraint. \code{llog}, \code{llogp} and
\code{llogpp} are used by both \code{el.test.wt3} and
\code{el.test.wt2}. In addition \code{logwelr}
is used by \code{el.test.wt3}.

\code{MMRtime} estimate the Mean and Median residual life time at age
               with right censored data.

\code{iter} is for perform one iteration of EM in the Buckley-James
            censored regression estimation.
\code{redistF} is for redistribution of probability, according to Fdist.
Used in BJtest().
}
\keyword{internal}
