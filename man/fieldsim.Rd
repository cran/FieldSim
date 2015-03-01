\encoding{latin1}
\name{fieldsim}
\alias{fieldsim}

\title{Simulate manifold indexed Gaussian field by the Fieldsim method}

\usage{fieldsim(process,Ne,nbNeighbor)}

\arguments{

  \item{process}{an S4 object process}
  \item{Ne}{a positive integer corresponding to the number of points to simulate with the accurate simulation step}
\item{nbNeighbor}{a positive integer (between 2 and 32) corresponding to the number of neighbors to use in the second refined step of the algorithm.}

}

\description{
 The function \code{fieldsim} yields simulation of sample path of a 
manifold indexed Gaussian field (or bridge) following the
procedure described in Brouste et al. (2007, 2010, 2014).
}

\value{The function returns in the slot \code{values} of the object process the values of the process on the manifold atlas}

\references{A. Brouste, J. Istas and S. Lambert-Lacroix (2007). On Gaussian random 
fields simulations. A. Brouste, J. Istas and S. Lambert-Lacroix (2010) On simulation of manifold indexed fractional Gaussian fields. A. Brouste, J. Istas and S. Lambert-Lacroix (2014) Fractional Gaussian bridges with the package \code{FieldSim}.}

\author{Alexandre Brouste (\url{http://perso.univ-lemans.fr/~abrouste/}) and
 Sophie Lambert-Lacroix (\url{http://membres-timc.imag.fr/Sophie.Lambert/}).}

\seealso{\code{\link{process-class}}, \code{\link{setProcess}}.}

\examples{
# Load FieldSim library
library(FieldSim)

# Fractional Brownian field on [0,1]^2
plane.fBm<-setProcess("fBm-plane",0.7)
str(plane.fBm)
fieldsim(plane.fBm)
plot(plane.fBm)

# Sphere indexed fractional Brownian field
#sphere.fBm<-setProcess("fBm-sphere",0.3)
#fieldsim(sphere.fBm)
#plot(sphere.fBm)

# Bridge associated to the Fractional Brownian field on [0,1]^2
#Gamma<-matrix(c(1,0,0,0,1,1,1,1,1,1/2,1/2,0.5),3,4)
#bridge.plane.fBm<-setProcess("bridge-fBm-plane",list(Gamma=Gamma,par=0.9))
#fieldsim(bridge.plane.fBm)
#plot(bridge.plane.fBm)

# Other examples can be found in the setProcess documentation.
}