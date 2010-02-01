\encoding{latin1}
\name{constructcovf}
\alias{constructcovf}
\title{Construct usual fractional covariance functions}
\usage{constructcovf(manifold,typeproc,H,F)}
\arguments{
\item{manifold}{a manifold;}
\item{typeproc}{the type of covariance, possible choice \code{"fBm"} or \code{"mBm"};}
\item{H}{Hurst parameter for fBm;}
\item{F}{Hurst function for mBm.}}
\description{The function \code{constructcovf} constructs the 
covariance function of usual fractional processes (fBm, mBm).}
\references{
A. Brouste, J. Istas and S. Lambert-Lacroix (2007). On Gaussian random 
fields simulations.
A. Brouste, J. Istas and S. Lambert-Lacroix (2010). On simulation of manifold indexed fractional Gaussian fields.}

\seealso{\code{\link{fieldsim}}.}
\author{
  Alexandre Brouste (\url{http://subaru.univ-lemans.fr/sciences/statist/pages_persos/Brouste/})
   and
  Sophie Lambert-Lacroix (\url{http://ljk.imag.fr/membres/Sophie.Lambert}).
}

\examples{
# Load FieldSim library
library(FieldSim)

plane<-setManifold("plane")

#Two dimensional fractional Brownian motion
R<-constructcovf(plane,"fBm",H=0.6)

#Two dimensional multifractional Brownian motion
F.2<-function(x){return(1/4+1/2*x[1])}
R.2<-constructcovf(plane,"mBm",F=F.2)

#For spherical and hyperboloid examples, see fieldsim documentation.
}
