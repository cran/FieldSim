\encoding{latin1}
\name{quadvar}
\alias{quadvar}

\title{Estimate the Hurst parameter of a plane indexed fractional Brownian field by the quadratic variations method}
\usage{
quadvar(process,parameter)
}

\arguments{
  \item{process}{a S4 object process;}
  \item{parameter}{parameter (in progress).}
}

\description{
 The function \code{quadvar} yields the estimation of the Hurst parameter of a fractional Brownian field
by the quadratic variations method in the plane case.
}

\details{
 The Hurst parameter of the fractal Brownian field is estimated by the procedure described
in Istas and Lang (1997).
}


\value{\item{H}{
a real in \eqn{]0,1[} that represents the estimate of the Hurst parameter of the fractional Brownian field.}
  }

\references{J. Istas and G. Lang (1997). Quadratic variations and estimation of the local
Holder index of a Gaussian process. Annales Institut Henri Poincare, 33,407-436.
}

\author{Alexandre Brouste (\url{http://perso.univ-lemans.fr/~abrouste/}) and
 Sophie Lambert-Lacroix (\url{http://membres-timc.imag.fr/Sophie.Lambert/}).}

\seealso{\code{\link{fieldsim}}, \code{\link{setProcess}}, \code{\link{setValues}}.}

\examples{
# load FieldSim library
library(FieldSim)

# Simulated Fractional Brownian field on [0,1]^2
plane.fBm<-setProcess("fBm-plane",0.7)
fieldsim(plane.fBm)
quadvar(plane.fBm)

# Simulated Multifractional Brownian field on [0,1]^2
funcH<-function(xi){0.3+xi[1]*0.6}
plane.mBm<-setProcess("mBm-plane",funcH)
fieldsim(plane.mBm)
quadvar(plane.mBm,parameter=list(point=c(0.5,0.5),h=0.2))
}
