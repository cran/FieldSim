\encoding{latin1}
\name{plot}
\alias{plot}
\title{Generic plotting of specific manifold indexed fractional Gaussian processes}
\usage{plot(x,y,...)}

\arguments{
  \item{x}{an object of class \code{process};}
  \item{y}{the type of the plot, possible choices are \code{"default"}, \code{"cloud"} or \code{"sun"}.}
   \item{...}{Other plot arguments}
}

\description{The function plots some of usual manifold indexed fractional Gaussian processes.}

\references{A. Brouste, J. Istas and S. Lambert-Lacroix (2010). On simulation of manifold indexed fractional Gaussian fields.}

\author{Alexandre Brouste (\url{http://perso.univ-lemans.fr/~abrouste/}) and
 Sophie Lambert-Lacroix (\url{http://membres-timc.imag.fr/Sophie.Lambert/}).}

\seealso{\code{\link{fieldsim}}, \code{\link{setProcess}}.}

\examples{
# Load FieldSim library
library(FieldSim)

# Fractional Brownian field on [0,1]^2
plane.fBm<-setProcess("fBm-plane",0.7)
str(plane.fBm)
fieldsim(plane.fBm)
plot(plane.fBm)

#The "cloud" plotting
plot(plane.fBm,"cloud")

#The "sun" plotting
plot(plane.fBm,"sun")

# Sphere indexed fractional Brownian field
#sphere.fBm<-setProcess("fBm-sphere",0.3)
#fieldsim(sphere.fBm)
#plot(sphere.fBm)
}


