\encoding{latin1}
\name{fieldsim}
\alias{fieldsim}

\title{Random field simulation by the FieldSim method}
\usage{fieldsim(S,Elevel=1,Rlevel=5,nbNeighbor=4)}

\arguments{
  \item{S}{a covariance function (defined on \eqn{[0,1]^4}) 
  of a Random field to simulate.}
  \item{Elevel}{a positive integer. \eqn{\code{Elevel}} is the
level associated with the regular space discretization for the first 
step of the algorithm (Accurate simulation step). The 
regular space discretization is of the following form: 
\eqn{[[0:2^{\code{Elevel}}]/2^{\code{Elevel}}]^2}.}
  \item{Rlevel}{a positive integer. \eqn{\code{Elevel}+\code{Rlevel}} is the
level associated with the regular space discretization for the second 
step of the algorithm (Refined simulation step). The 
regular space discretization is of the following form: 
\eqn{[[0:2^{\code{Elevel}+\code{Rlevel}}]/2^{\code{Elevel}+\code{Rlevel}}]^2}.}
\item{nbNeighbor}{a positive integer. \code{nbNeighbor} must be between 1 and 32.
\code{nbNeighbor} is the number of neighbors to use in the second step 
of the algorithm.}
}
\description{
 The function \code{fieldsim} yields discretization of sample path of a Gaussian field following the
procedure described in Brouste et al. (2007). 
  }

\details{The function \code{fieldsim} yields discretization of sample path of a Gaussian field
associated with the covariance function given by \code{R}.
 The subspace [0,1]x[0,1] is discretized in a regular space discretization of size
\eqn{(2^{\code{Elevel}+\code{Rlevel}}+1)^2}. At each point of the grid, the Gaussian field is simulated using the
procedure described in Brouste et al. (2007).
}



\value{
  A list with the following components:
  \item{Zrow}{the vector of length \eqn{2^{\code{Elevel}+\code{Rlevel}}+1} containing the discretization 
  of the x axis.} 
  \item{Zcol}{the vector of length \eqn{2^{\code{Elevel}+\code{Rlevel}}+1} containing the discretization 
  of the y axis.}
  \item{Z}{the matrix of size \eqn{(2^{\code{Elevel}+\code{Rlevel}}+1)}x\eqn{(2^{\code{Elevel}+\code{Rlevel}}+1)}
  in such a way \code{Z}[i,j] containing the value of the process at point (\code{Zrow}[i],\code{Zcol}[j])}
  \item{time}{the CPU time}
  }

\references{A. Brouste, J. Istas and S. Lambert-Lacroix (2007). On Gaussian random 
fields simulations.}

\author{Alexandre Brouste (\url{http://ljk.imag.fr/membres/Alexandre.Brouste}) and
 Sophie Lambert-Lacroix (\url{http://ljk.imag.fr/membres/Sophie.Lambert}).}

\seealso{\code{\link{quadvar}}.}

\examples{
# load FieldSim library
library(FieldSim)

#Example 1: Fractional Brownian Field
R<-function(x,H=0.9){1/2*((x[1]^2+x[2]^2)^(H)+(x[3]^2+x[4]^2)^(H)-((x[1]-x[3])^2+(x[2]-x[4])^2)^(H))}
res<- fieldsim(R,Elevel=1,Rlevel=5,nbNeighbor=4)

# Plot 
x <- res$Zrow
y <- res$Zcol
z <- res$Z
persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")

#Example 1: Multifractional Brownian Field

F<-function(y){0.4*y+0.5}

R<-function(x,Fun=F){
H1<-Fun(x[1])
H2<-Fun(x[3])
alpha<-1/2*(H1+H2)
C2D(alpha)^2/(2*C2D(H1)*C2D(H2))*((x[1]^2+x[2]^2)^(alpha)+(x[3]^2+x[4]^2)^(alpha)-((x[1]-x[3])^2+(x[2]-x[4])^2)^(alpha))
}

res<- fieldsim(R,Elevel=1,Rlevel=5,nbNeighbor=4)
# Plot 
x <- res$Zrow
y <- res$Zcol
z <- res$Z
persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")
}
\keyword{}
