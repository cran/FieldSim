\encoding{latin1}
\name{setManifold}
\alias{setManifold}

\title{Set a S4 manifold object}

\usage{setManifold(name,atlas,gridtype,distance,origin)}

\arguments{
  \item{name}{name of the manifold (type character);}
  \item{atlas}{atlas of the manifold (type matrix);}
  \item{gridtype}{is the  grid type (a character string) to plotting;}
  \item{distance}{distance on the manifold (type function);}
  \item{origin}{origin of the manifold (type matrix).}
}

\description{The function sets an object of class manifold.}
\value{ An object of class manifold with the 5 slots name, atlas, gridtype, distance and origin.}

\references{A. Brouste, J. Istas and S. Lambert-Lacroix (2010). On simulation of manifold indexed fractional Gaussian fields.}

\author{Alexandre Brouste (\url{http://perso.univ-lemans.fr/~abrouste/}) and Sophie Lambert-Lacroix (\url{http://membres-timc.imag.fr/Sophie.Lambert/}).}

\seealso{\code{\link{setProcess}}.}


\examples{
# Load FieldSim library
library(FieldSim)

# Example 1: User manifold
name1<-"plane1"
mesh<-seq(from=0,to=1,length=16)
atlas1<-rbind(rep(mesh,each=16),rep(mesh,16))
d1<-function(xi,xj){return(sqrt(t(xi-xj)\%*\%(xi-xj)))}
origin1<-rbind(0,0)
manifold1<-setManifold(name=name1, atlas=atlas1, distance=d1, origin=origin1)
str(manifold1)

#Example 2: The "line" manifold
line<-setManifold("line")
str(line)

#Example 3: The "plane" manifold
plane<-setManifold("plane")
str(plane)

#Example 4: The "sphere" manifold
sphere<-setManifold("sphere")
str(sphere)

#Example 5: The "hyperboloid" manifold
hyper<-setManifold("hyperboloid")
str(hyper)
}


