\name{emplik-package}
\alias{emplik-package}
\alias{emplik}
\docType{package}
\title{
Empirical likelihood for mean functional/hazard functional with possibly censored data. 
}
\description{
Empirical likelihood ratio tests and confidence intervals for means/hazards
 	from possibly censored and/or truncated data. In particular, empirical likelihood for the
	Kaplan-Meier/Nelson-Aalen estimators. Now does AFT regression too.
}
\details{
For non-censored data and mean parameters, use \code{el.test( )}.

For censored data and mean parameters, use \code{el.cen.EM2( )}.

For censored data and hazard parameters, use \code{emplikH1.test( )} [Poisson type likelihood];
     use \code{emplikH.disc( )} [binomial type likelihood].

For constructing confidence intervals from likelihood ratio, use \code{findUL( )}.
% An overview of how to use the package, including the most important ~~
% functions ~~
}
\author{
Mai Zhou. (el.test( ) is adapted from Owen's splus code elm. Yifan Yang for some C code.)

Maintainer: Mai Zhou <maizhou@gmail.com>
}
\references{
Zhou, M. (2016).
Empirical Likelihood Method in Survival Analysis.
(Chapman and Hall/CRC Biostatistics Series)
CRC press 2016
}
\keyword{package}
\seealso{Another R package kmc for possible faster results when testing mean functional with right censored data.}