\encoding{latin1}
\name{midpoint}
\alias{midpoint}

\title{Fractional Brownian field simulation by the midpoint displacement method}
\usage{
midpoint(process)
}
\arguments{
  \item{process}{an object of class process (namely an FBm).}
  }

\description{
 The function \code{midpoint} yields simulation of sample path of a fractional Brownian field
by the midpoint displacement method.
  }

\details{
 The subspace [0,1] x [0,1] is discretized in a regular space discretization of size
\eqn{(2^{\code{nblevel}}+1)^2}. At each point of the grid, the fractional Brownian
field is simulated using the midpoint displacement method described for example
in Fournier et al. (1982).
}

\value{
 an object of class process with the simulated sample path in the corresponding slot \code{values}.
  }

\references{A. Fournier, D. Fussel and L. Carpenter
(1982) Computer rendering of stochastic model, Communication of the AMC, 25,
371-384.

H.O. Peitgen and D. Saupe (1998) The science of fractal images, Springer-Verlag.

R.F. Voss (1985) Random fractal forgeries. NATO ASI Series, F17, 805-835.
}

\author{Alexandre Brouste (\url{http://perso.univ-lemans.fr/~abrouste/}) and
 Sophie Lambert-Lacroix (\url{http://membres-timc.imag.fr/Sophie.Lambert/}).}

\seealso{\code{\link{fieldsim}}.}

\examples{
# load FieldSim library
library(FieldSim)

plane.fBm<-setProcess("fBm-plane",0.9)
midpoint(plane.fBm)
plot(plane.fBm)
}
