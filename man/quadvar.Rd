\encoding{latin1}
\name{quadvar}
\alias{quadvar}

\title{Estimation of the Hurst parameter of a fractional Brownian field
            by the quadratic variations method}
\usage{
quadvar(Z)
}
\arguments{
  \item{Z}{a matrix of size (n x n). \code{Z} is the matrix associated with the sample path of one
fractal Brownian field. n must be of the form \eqn{2^{J}+1} where \eqn{J} is a positive integer.
So \code{Z}[i,j] is the value of the process at the point \eqn{((i-1)2^{-J},(j-1)2^{-J})}.
  }
}
\description{
 The function \code{quadvar} yields the estimation of the Hurst parameter of a fractional Brownian field
by the quadratic variations method.
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
  Alexandre Brouste (\url{http://ljk.imag.fr/membres/Alexandre.Brouste})
   and
  Sophie Lambert-Lacroix (\url{http://ljk.imag.fr/membres/Sophie.Lambert}).
}

\seealso{\code{\link{fieldsim}}, \code{\link{midpoint}}.}

\examples{
# load FieldSim library
library(FieldSim)

# Simulation
# H=0.3
res <- midpoint(H=0.3,nblevel=8)
# Estimation
quadvar(Z=res$Z)

}
\keyword{}
