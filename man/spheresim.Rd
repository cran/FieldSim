\encoding{latin1}
\name{spheresim}
\alias{spheresim}

\title{Random spherical field simulation by the mifieldsim method on visual grid}
\usage{spheresim(R,Ne=100,Nr=100,Ng=100,nbNeighbor=4)}

\arguments{
  \item{R}{a covariance function (defined on the sphere) 
  of a Random spherical field to simulate.}
  \item{Ne}{a positive integer. \eqn{\code{Ne}} is the number 
     of simulation points associated with the uniform distributed discretization 
     of the sphere for the first step of the algorithm (Accurate simulation step)}
  \item{Nr}{a positive integer. \eqn{\code{Nr}} is the number of simulation points 
  associated with the uniform distributed discretization 
   of the sphere for the second step of the algorithm (Refined simulation step).}
  \item{Ng}{a positive integer. \eqn{\code{Nr}} is the number of simulation points 
    associated with the visual grid discretization 
    of the sphere for the third step of the algorithm 
    (Visual refined simulation step).}
 \item{nbNeighbor}{a positive integer. \code{nbNeighbor} must be between 1 and 32.
   \code{nbNeighbor} is the number of neighbors to use in the second step 
    of the algorithm.}
}
\description{
 The function \code{spheresim} yields discretization of sample path of a 
 Gaussian spherical field following the
 procedure described in Brouste et al. (2009). 
  }

\details{The function \code{spheresim} yields discretization of sample path of a 
Gaussian spherical field associated with the covariance function given by \code{R}.}

\value{
  A list with the following components:
  \item{X}{the vector of length at more \eqn{\code{Ne}+\code{Nr}+6\code{Nr}^2} 
  containing the discretization of the x axis.} 
  \item{Y}{the vector of length at more \eqn{\code{Ne}+\code{Nr}+6\code{Nr}^2} 
  containing the discretization of the y axis.}
  \item{Z}{the vector of length at more \eqn{\code{Ne}+\code{Nr}+6\code{Nr}^2} 
  containing the discretization of the z axis.}
  \item{W}{the vector of length at more \eqn{\code{Ne}+\code{Nr}+6\code{Nr}^2} 
  containing the value of the simulated field at point 
  (\code{X}[n],\code{Y}[n],\code{Z}[n])}
  \item{W1,...,W6}{the matrices of of size \eqn{\code{Ng}^2} 
   that give values of the simulated spherical field at the points 
   of the visual grid}
  \item{time}{the CPU time}
  }

\references{A. Brouste, J. Istas and S. Lambert-Lacroix (2009). On simulation of manifold indexed fractional Gaussian fields.}

\author{Alexandre Brouste (\url{http://subaru.univ-lemans.fr/sciences/statist/pages_persos/Brouste/}) and
 Sophie Lambert-Lacroix (\url{http://ljk.imag.fr/membres/Sophie.Lambert}).}

\seealso{\code{\link{spheresimgrid}},\code{\link{fieldsim}},\code{\link{hypersim}}.}

\examples{
# load FieldSim library
library(FieldSim)

d<-function(x){ #Distance on the sphere
u <- x[1]*x[4]+x[2]*x[5]+x[3]*x[6]
if (u<(-1))
u<--1
if (u>1)
u<-1
acos(u)
}

# Example 1 : Fractional spherical field with RS1 covariance function

RS1<-function(x){
H<-0.45	# H can vary from 0 to 0.5
1/2*(d(c(1,0,0,x[1:3]))^{2*H}+d(c(1,0,0,x[4:6]))^{2*H}-d(x)^{2*H})
}

res1<- spheresim(RS1,Ne=100,Nr=1000,Ng=25,nbNeighbor=4)

library(rgl)
library(RColorBrewer)
printsphere(res1)

# Example 2 : Fractional spherical field with RS2 covariance function

RS2<-function(x){
H<-0.45	# H can vary from 0 to 0.5
exp(-d(x)^{2*H})
}

res2<- spheresim(RS2,Ne=100,Nr=1000,Ng=25,nbNeighbor=4)
printsphere(res2)
}
