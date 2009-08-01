\encoding{latin1}
\name{midpoint}
\alias{midpoint}

\title{Fractional Brownian field simulation by the midpoint displacement method}
\usage{
midpoint(H,nblevel=8)
}
\arguments{
  \item{H}{a real in ]0,1[. \code{H} is the Hurst parameter of the fractional Brownian field to simulate.}
  \item{nblevel}{a positive integer. \eqn{\code{nblevel}} is the
level associated with the regular space discretization of the following form:
\eqn{[[0:2^{\code{nblevel}}]/2^{\code{nblevel}}]^2}.}
  }
\description{
 The function \code{midpoint} yields discretization of sample path of a fractional Brownian field
by the midpoint displacement method.
  }

\details{
 The subspace [0,1] x [0,1] is discretized in a regular space discretization of size
\eqn{(2^{\code{nblevel}}+1)^2}. At each point of the grid, the fractional Brownian
field is simulated using the midpoint displacement method described for example
in Fournier et al. (1982).
}



\value{
  A list with the following components:
  \item{Zrow}{the vector of length \eqn{2^{\code{nblevel}}+1} containing the discretization
  of the x axis.}
  \item{Zcol}{the vector of length \eqn{2^{\code{nblevel}}+1} containing the discretization
  of the y axis.}
  \item{Z}{the matrix of size \eqn{(2^{\code{nblevel}}+1)}x\eqn{(2^{\code{nblevel}}+1)}
  in such a way \code{Z}[i,j] containing the value of the process at point (\code{Zrow}[i],\code{Zcol}[j])}
  \item{time}{the CPU time}
  }

\references{A. Fournier, D. Fussel and L. Carpenter
(1982). Computer rendering of stochastic models. Communication of the AMC *25*,
371-384.

H.O. Peitgen and D. Saupe (1998). The science of fractal images. Springer Verlag.

R.F. Voss (1985). Random fractal forgeries. NATO ASI Series *F17*, 805-835.
}

\author{
  Alexandre Brouste (\url{http://subaru.univ-lemans.fr/sciences/statist/pages_persos/Brouste/})
   and
  Sophie Lambert-Lacroix (\url{http://ljk.imag.fr/membres/Sophie.Lambert}).
}

\seealso{\code{\link{quadvar}},\code{\link{fieldsim}}}

\examples{
# load FieldSim library
library(FieldSim)

# Simulation
# H=0.3
res <- midpoint(H=0.3,nblevel=8)
# Plot
x <- res$Zrow
y <- res$Zcol
z <- res$Z
persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")
}
