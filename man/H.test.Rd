\encoding{latin1}
\name{H.test}
\alias{H.test}

\title{Performs Hurst parameter test for random Gaussian field}
\usage{H.test(Z,H,alternative="two.sided",conf.level=0.95)}

\arguments{
  \item{Z}{a matrix of size (n x n). \code{Z} is the matrix associated with the sample path of one
fractional Brownian field. n must be of the form \eqn{2^{J}+1} where \eqn{J} is a positive integer.
So \code{Z}[i,j] is the value of the process at the point \eqn{((i-1)2^{-J},(j-1)2^{-J})}.}
  \item{H}{a real in ]0,1[. \code{H} is the Hurst parameter of the fractional Brownian 
  field to be tested.}
  \item{alternative}{a character string specifying the alternative hypothesis,
          must be one of '"two.sided"' (default), '"greater"' or
          '"less"'.}
  \item{conf.level}{confidence level of the interval}
  }
}

\description{ The function \code{H.test} performs the test "the Hurst parameter of Z is 
equal to H" against the alternative hypothesis.}
\details{The test is build using the asymptotic distribution of the quadratic variations 
given in Brouste et al. (2006)}

\value{
  A list with the following components:
  \item{statistic}{the value of the test statistic,}
  \item{p.value}{the p-value for the test,}
  \item{conf.int}{a confidence interval for the Hurst parameter,}
  \item{gamma2}{the constant of the variance used in the test.}
}
\references{A. Brouste, J. Istas and S. Lambert-Lacroix (2007). 
   On Gaussian random fields simulations.
}

\author{Alexandre Brouste (\url{http://ljk.imag.fr/membres/Alexandre.Brouste}) 
and Sophie Lambert-Lacroix (\url{http://ljk.imag.fr/membres/Sophie.Lambert}).}

\seealso{\code{\link{fieldsim}}, \code{\link{midpoint}}.}

\examples{
# load FieldSim library
library(FieldSim)

# Simulation
# H=0.3
res <- midpoint(H=0.3,nblevel=8)
# Test "H=0.3" against "H!=3"
H.test(res$Z,H=0.3)
}
\keyword{}
