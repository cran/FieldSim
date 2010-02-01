\encoding{latin1}
\name{quadvar}
\alias{quadvar}

\title{Estimate the Hurst parameter of a plane indexed fractional Brownian field by the quadratic variations method}
\usage{
quadvar(manifold,res)
}

\arguments{
  \item{manifold}{a S4 object manifold.}
  \item{res}{sample path of the field on the manifold atlas}
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
Holder index of a Gaussian process. Annales Institut Henri Poincaré *33* 407-436.
}

\author{
  Alexandre Brouste (\url{http://subaru.univ-lemans.fr/sciences/statist/pages_persos/Brouste/})
   and
  Sophie Lambert-Lacroix (\url{http://ljk.imag.fr/membres/Sophie.Lambert}).
}

\seealso{\code{\link{fieldsim}}.}

\examples{
# load FieldSim library
library(FieldSim)

plane<-setManifold("plane")
R<-constructcovf(plane,"fBm", H=0.4)
res<-fieldsim(plane,R,Ne=80,nbNeighbor=15)

quadvar(plane,res)

}
