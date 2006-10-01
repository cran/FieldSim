\encoding{latin1}
\name{locquadvar}
\alias{locquadvar}

\title{Estimation of the multifractional function of the multifractalional Brownian field
            by the localized quadratic variations method}

\usage{
locquadvar(Z,t,h)
}
\arguments{
  \item{Z}{a matrix of size (n x n). \code{Z} is the matrix associated with the sample path of one
fractal Brownian field. n must be of the form \eqn{2^{J}+1} where \eqn{J} is a positive integer.
So \code{Z}[i,j] is the value of the process at the point \eqn{((i-1)2^{-J},(j-1)2^{-J})}.}
  \item{t}{a vector of length 2. \code{t} is the point of the square [0,1] x [0,1]
  where multifractional function must be estimated.}
  \item{h}{a positive real. \code{h} defines the neighborhood used to estimate multifractional
  function via the following form \eqn{((i-1)2^{-J},(j-1)2^{-J})}, \eqn{i,j=1,...,n} such that
  \eqn{abs(2^{-iJ}-t[1])<=h} and \eqn{abs(2^{-jJ}-t[2])<=h}.}
}
\description{
 The function \code{locquadvar} yields the estimation of the multifractional function at point \code{t}
 of the multifractional Brownian field by the localized quadratic variations method.
  }

\details{
 The multifractional function at point \code{t} of the multi-fractal Brownian field
 is estimated by the procedure described in Lacaux (2004).
}


\value{
  A list with the following components:
  \item{H}{
a real in \eqn{]0,1[} that represents the estimate of
the multifractional function at point \code{t}.}
  }

\references{C. Lacaux (2004). Real Harmonizable Multifractal Levy Motions.
Annales Institut Henri Poincaré, Prob. Stat., *40* 259-277.
}

\author{Alexandre Brouste (\url{http://ljk.imag.fr/membres/Alexandre.Brouste}) and
  Sophie Lambert-Lacroix (\url{http://ljk.imag.fr/membres/Sophie.Lambert}).}

\seealso{\code{\link{fieldsim}}.}

\examples{
# load FieldSim library
library(FieldSim)

#Field simulation
F<-function(x){
H1<-(0.4*x[1]+0.5)
H2<-(0.4*x[3]+0.5)
out<-list(H1=H1,H2=H2,H=(H1+H2)/2)
return(out)
}

R<-function(x){
a<-F(x)
alpha<-a$H
alpha1<-a$H1
alpha2<-a$H2

C2D(alpha)^2/(2*C2D(alpha1)*C2D(alpha2))*((x[1]^2+x[2]^2)^(alpha)+(x[3]^2+x[4]^2)^(alpha)-((x[1]-x[3])^2+(x[2]-x[4])^2)^(alpha))
}

res <- fieldsim(R,Elevel=1,Rlevel=5,nbNeighbor=4)

# Estimation
locquadvar(Z=res$Z,t=c(0.5,0.5),h=0.125)
}
\keyword{}
