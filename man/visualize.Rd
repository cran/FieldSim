\encoding{latin1}
\name{visualize}
\alias{visualize}
\title{Visualize some of specific manifold indexed fractional Gaussian process}
\usage{visualize(manifold,res,typeplot="default")}
\arguments{
  \item{manifold}{a manifold;}
  \item{res}{simulation heights of the process on the atlas of the manifold (type vector);}
  \item{typeplot}{the type of the plot, possible choice \code{"default"}, \code{"cloud"} or \code{"sun"};;}}
\description{The function plots some of usual manifold indexed fractional Gaussian processes.}

\references{A. Brouste, J. Istas and S. Lambert-Lacroix (2010). On simulation of manifold indexed fractional Gaussian fields.}

\author{Alexandre Brouste (\url{http://subaru.univ-lemans.fr/sciences/statist/pages_persos/Brouste/}) and
 Sophie Lambert-Lacroix (\url{http://ljk.imag.fr/membres/Sophie.Lambert}).}

\seealso{\code{\link{setManifold}}, \code{\link{fieldsim}}}

\examples{
# Load FieldSim library
library(FieldSim)

#Example 1
plane<-setManifold("plane")
R<-constructcovf(plane,"fBm",H=0.6)
res<-fieldsim(plane,R,50)

visualize(plane,res)

#Example 2: The "cloud" plotting
visualize(plane,res,"cloud")

#Example 3: The "sun" plotting
visualize(plane,res,"sun")
}


