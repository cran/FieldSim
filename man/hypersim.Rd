\encoding{latin1}
\name{hypersim}
\alias{hypersim}

\title{Random hyperboloïd indexed field simulation by the mifieldsim method on visual grid}
\usage{hypersim(R,Ne=100,Nr=100,Ng=100,nbNeighbor=4)}

\arguments{
  \item{R}{a covariance function (defined on the hyperboloïd) of a Random 
		hyperboloïd indexed field to simulate.}
  \item{Ne}{a positive integer. \eqn{\code{Ne}} is the number 
     	of simulation points associated with the uniform distributed discretization 
     	of the hyperboloïd for the first step of the algorithm (Accurate simulation step)}
  \item{Nr}{a positive integer. \eqn{\code{Nr}} is the number of simulation points 
  associated with the uniform distributed discretization 
   of the sphere for the second step of the algorithm (Refined simulation step).}
  \item{Ng}{a positive integer. \eqn{\code{Nr}} is the number of simulation points 
    associated with the visual grid discretization 
    of the hyperboloïd for the third step of the algorithm 
    (Visual refined simulation step).}
 \item{nbNeighbor}{a positive integer. \code{nbNeighbor} must be between 1 and 32.
   \code{nbNeighbor} is the number of neighbors to use in the second step 
    of the algorithm.}
}

\description{
 The function \code{hypersim} yields discretization of sample path of a 
 Gaussian hyperboloïd indexed field following the procedure described in 
Brouste et al. (2009).
}

\details{The function \code{hypersim} yields discretization of sample path 
of a Gaussian hyperboloïd indexed field associated with the covariance function given 
by \code{R}.}

\value{
  A list with the following components:
  \item{X}{the vector of length at more \eqn{\code{Ne}+\code{Nr}+\code{Nr}^2} 
  containing the discretization of the x axis.} 
  \item{Y}{the vector of length at more \eqn{\code{Ne}+\code{Nr}+\code{Nr}^2} 
  containing the discretization of the y axis.}
  \item{Z}{the vector of length at more \eqn{\code{Ne}+\code{Nr}+\code{Nr}^2} 
  containing the discretization of the z axis.}
  \item{W}{the vector of length at more \eqn{\code{Ne}+\code{Nr}+\code{Nr}^2} 
  containing the value of the simulated field at point 
  (\code{X}[n],\code{Y}[n],\code{Z}[n])}
  \item{W1}{the matrice of size \eqn{\code{Ng}^2} 
   that give values of the simulated hyperboloïd indexed field at the points 
   of the visual grid}
  \item{time}{the CPU time}
  }

\references{A. Brouste, J. Istas and S. Lambert-Lacroix (2009). On simulation of manifold indexed fractional Gaussian fields.}

\author{Alexandre Brouste (\url{http://subaru.univ-lemans.fr/sciences/statist/pages_persos/Brouste/}) and
 Sophie Lambert-Lacroix (\url{http://ljk.imag.fr/membres/Sophie.Lambert}).}

\seealso{\code{\link{hypersimgrid}},\code{\link{fieldsim}},\code{\link{spheresim}}.}

\examples{
# Load FieldSim library
library(FieldSim)

d<-function(x){    #Distance on the hyperboloïd
u <- -x[1]*x[4]-x[2]*x[5]+x[3]*x[6]
if (u<1){u<-1}
acosh(u)
}


#Example 1: Hyperboloïd indexed Brownian fractional Field with RH1 covariance function

RH1<-function(x){
H<-0.45          # H can vary from 0 to 0.5
1/2*(d(c(0,0,1,x[1:3]))^{2*H}+d(c(0,0,1,x[4:6]))^{2*H}-d(x)^{2*H})
}

resh1<- hypersim(RH1,Ne=100,Nr=1000,Ng=50,nbNeighbor=4)

library(rgl)
library(RColorBrewer)
printhyper(resh1)

#Example 2: Hyperboloïd indexed Brownian Field with RH4 covariance function

RH4<-function(x){
H<-0.45
1/(1+d(x)^{2*H})
}

resh4<- hypersim(RH4,Ne=100,Nr=1000,Ng=50,nbNeighbor=4)
printhyper(resh4)
}
