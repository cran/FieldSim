\encoding{latin1}
\name{hypersimgrid}
\alias{hypersimgrid}

\title{Random hyperboloid indexed field simulation by the mifieldsim 
method on the user grid}

\usage{hypersimgrid(R,Ne=100,Nr=100,nbNeighbor=4,S)}

\arguments{
  \item{R}{a covariance function (defined on the sphere) 
  of a Random spherical field to simulate.}
  \item{Ne}{a positive integer. \eqn{\code{Ne}} is the number 
     of simulation points associated with the uniform distributed discretization 
     of the hyperboloïd for the first step of the algorithm (Accurate simulation step)}
  \item{Nr}{a positive integer. \eqn{\code{Nr}} is the number of simulation points 
  associated with the uniform distributed discretization 
   of the hyperboloïd for the second step of the algorithm (Refined simulation step).}
 \item{nbNeighbor}{a positive integer. \code{nbNeighbor} must be between 1 and 32.
   \code{nbNeighbor} is the number of neighbors to use in the second step 
    of the algorithm.}
  \item{S}{a matrix with 3 rows and \code{N} columns. (\code{S}[1,n],\code{S}[2,n],\code{S}[3,n]) 
 is the point of the hyperboloïd where the field must be simulated.}
}
\description{
 The function \code{hypersimgrid} yields discretization of sample path of a 
 Gaussian hyperboloïd indexed field on the \code{S} grid following the
 procedure described in Brouste et al. (2009). 
  }

\details{The function \code{hypersimgrid} yields discretization of sample path of a 
Gaussian hyperboloid indexed field associated with the covariance function given 
by \code{R}.}

\value{
  A list with the following components:
  \item{S}{a matrix with 3 rows and \code{N} columns. (\code{S}[1,n],\code{S}[2,n],\code{S}[3,n]) 
 is the point where the field has been simulated.}
  \item{Wg}{the vector of length at more \code{N} 
  containing the value of the simulated field at point 
  (\code{S}[1,n],\code{S}[2,n],\code{S}[3,n])}
  \item{time}{the CPU time}
  }

\references{A. Brouste, J. Istas and S. Lambert-Lacroix (2009). On simulation of manifold indexed fractional Gaussian fields.}

\author{Alexandre Brouste (\url{http://subaru.univ-lemans.fr/sciences/statist/pages_persos/Brouste/}) and
 Sophie Lambert-Lacroix (\url{http://ljk.imag.fr/membres/Sophie.Lambert/}).}

\seealso{\code{\link{hypersim}}}

\examples{
# Load FieldSim library
library(FieldSim)

## Define the grid
Z <- 2
n <- 30 
phi <- Z/n*(0:(n-1))
theta <- 2*pi/n*(0:(n-1))
S <- numeric(0)
for (kk in 1:(n-1)){
S <- cbind(S,rbind(cos(theta[kk])*sinh(phi),sin(theta[kk])*sinh(phi),cosh(phi)))
}

## Define the autocovariance function
d<-function(x){
u <- -x[1]*x[4]-x[2]*x[5]+x[3]*x[6]
if (u<1){u<-1}
acosh(u)
}
R<-function(x){
H<-0.3
1/2*(d(c(0,0,1,x[1:3]))^{2*H}+d(c(0,0,1,x[4:6]))^{2*H}-d(x)^{2*H})
}

Wg <- hypersimgrid(R,Ne=50,Nr=100,nbNeighbor=4,S)$Wg

}
