\encoding{latin1}
\name{C2D}
\alias{C2D}
\title{Constant for multifractional Brownian field simulation}
\usage{C2D(H)}

\arguments{\item{H}{a real in ]0,1[. \code{H} is the Hurst parameter of the 
constant to simulate.}}
\description{The function \code{C2D} yields the value of the constant in the covariance 
function of the multifractional Brownian field (for the exact expression see 
Brouste et al. (2007) Section 3.2).}

\value{The value of the constant}

\author{
  Alexandre Brouste (\url{http://ljk.imag.fr/membres/Alexandre.Brouste})
   and
  Sophie Lambert-Lacroix (\url{http://ljk.imag.fr/membres/Sophie.Lambert}).
}

\seealso{\code{\link{fieldsim}}.}
\references{A. Brouste, J. Istas and S. Lambert-Lacroix (2007). On Gaussian random 
fields simulations.}

\examples{
# load FieldSim library
library(FieldSim)

# Constant
H<-0.3
C2D(H)
}
\keyword{}
