\name{emplik-package}
\alias{emplik-package}
\alias{emplik}
\docType{package}
\title{
Empirical likelihood for mean functional/hazard functional with possibly censored data. 
}
\description{
Empirical likelihood ratio tests and confidence intervals for means/hazards
 	from possibly censored and/or truncated data. 
	Now does regression too. The package contains some C code.
}
\details{
\tabular{ll}{
Package: \tab emplik      \cr
Type:    \tab Package     \cr
Version: \tab 1.1-0       \cr
Date:    \tab 2020-05-18  \cr
License: \tab GPL (>= 2)  \cr
}
For non-censored data and mean parameters, use \code{el.test( )}.

For censored data and mean parameters, use \code{el.cen.EM2( )}.

For censored data and hazard parameter, use \code{emplikH1.test( )} [Poisson type likelihood];
     use \code{emplikH.disc( )} [binomial type likelihood].

For constructing confidence interval, use \code{findUL( )}.
% An overview of how to use the package, including the most important ~~
% functions ~~
}
\author{
Mai Zhou (el.test is adapted from Owen's splus code. Yifan Yang for some C code.)

Maintainer: Mai Zhou <mai@ms.uky.edu> <maizhou@gmail.com>
}
\references{
Zhou, M. (2016).
Empirical Likelihood Method in Survival Analysis.
(Chapman and Hall/CRC Biostatistics Series)
CRC press 2016
}
\keyword{package}
\seealso{Another R package kmc for possible faster result for testing of mean with right censored data.}