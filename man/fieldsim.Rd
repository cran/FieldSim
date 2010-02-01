\encoding{latin1}
\name{fieldsim}
\alias{fieldsim}

\title{Simulate manifold indexed Gaussian field by the fieldsim method}

\usage{fieldsim(manifold,R,Ne,nbNeighbor)}

\arguments{
  \item{manifold}{an S4 object manifold}
  \item{R}{a covariance function of the Gaussian random field to simulate}
  \item{Ne}{a positive integer corresponding to the number of points to simulate with the accurate simulation step}
\item{nbNeighbor}{a positive integer (between 1 and 32) corresponding to the number of neighbors to use in the second refined step 
of the algorithm.}
}

\description{
 The function \code{fieldsim} yields discretization of sample path of a 
manifold indexed Gaussian field following the
procedure described in Brouste et al. (2007, 2010).
}

\value{The function returns the vector of the values of the process on the manifold atlas}

\references{A. Brouste, J. Istas and S. Lambert-Lacroix (2007). On Gaussian random 
fields simulations. A. Brouste, J. Istas and S. Lambert-Lacroix (2010) On simulation of manifold indexed fractional Gaussian fields. }

\author{Alexandre Brouste (\url{http://subaru.univ-lemans.fr/sciences/statist/pages_persos/Brouste/}) and
 Sophie Lambert-Lacroix (\url{http://ljk.imag.fr/membres/Sophie.Lambert}).}

\seealso{\code{\link{quadvar}},\code{\link{manifold-class}}}

\examples{
# Load FieldSim library
library(FieldSim)

# Sphere indexed fractional Brownian field

sphere<-setManifold("sphere")
R.S.1 <- constructcovf(sphere, "fBm", H = 0.4)

S.u <- constructgrid(sphere, "random", 10)
S.g <- constructgrid(sphere, "visualization", 12)
simulationgrid<-cbind(S.u,S.g)
sphere@atlas <- simulationgrid

resS <- fieldsim(sphere, R.S.1, Ne = 80, nbNeighbor = 15)

sphere@atlas<-S.g
res<-resS[(dim(S.u)[2]+1):length(resS)]
visualize(sphere,resS)
}