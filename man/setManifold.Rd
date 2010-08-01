\encoding{latin1}
\name{setManifold}
\alias{setManifold}
\title{Set a S4 manifold object}
\usage{setManifold(name,atlas,distance,origin)}
\arguments{
  \item{name}{name of the manifold (type character);}
  \item{atlas}{atlas of the manifold (type matrix);}
  \item{distance}{distance on the manifold (type function);}
  \item{origin}{origin of the manifold (type matrix).}}
\description{The function sets an object of class manifold.}
\value{ An object of class manifold with the 4 slots name, atlas, distance and origin.}

\references{
A. Brouste, J. Istas and S. Lambert-Lacroix (2010). On Simulation of Manifold Indexed Fractional Gaussian Fields. Journal of Statistical Software, 36(4), 1--14. URL \url{http://www.jstatsoft.org/v36/i04/}.
}

\author{Alexandre Brouste (\url{http://subaru.univ-lemans.fr/sciences/statist/pages_persos/Brouste/}) and
 Sophie Lambert-Lacroix (\url{http://ljk.imag.fr/membres/Sophie.Lambert}).}

\seealso{\code{\link{constructgrid}}, \code{\link{constructcovf}}, \code{\link{fieldsim}}, \code{\link{visualize}}}


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

#Example 2: The "plane" manifold
plane<-setManifold("plane")
str(plane)

#Example 3: The "sphere" manifold
sphere<-setManifold("sphere")
str(sphere)

#Example 4: The "hyperboloid" manifold
hyper<-setManifold("hyperboloid")
str(hyper)
}


