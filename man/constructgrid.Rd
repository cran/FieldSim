\encoding{latin1}
\name{constructgrid}
\alias{constructgrid}
\title{Construct usual grids on some specific manifolds}
\usage{constructgrid(manifold,typegrid,Ng)}
\arguments{
\item{manifold}{a manifold;}
\item{typegrid}{the type of the grid, possible choice \code{"regular"}, \code{"random"} or \code{"visualization"};}
\item{Ng}{parameter of the size of the grid, see details.}}

\description{The function \code{constructgrid} constructs usual grids on manifold.}

\references{
A. Brouste, J. Istas and S. Lambert-Lacroix (2010). On Simulation of Manifold Indexed Fractional Gaussian Fields. Journal of Statistical Software, 36(4), 1--14. URL \url{http://www.jstatsoft.org/v36/i04/}.
}


\details{ 
We list here the different implemented grids. For 
\code{manifold@name=="plane"} we have the \code{typegrid=="regular"} grid (with the parameter \code{Ng} returns a regular grid on [0,1]x[0,1] of size \code{Ng*Ng}), the \code{typegrid=="random"} grid (uniform random choice of the both coordinates on [0,1], grid of size \code{Ng*Ng}) and the \code{typegrid=="random"|typegrid=="visualization"} grid, of size \code{ (2^Ng +1)*(2^Ng +1)} composed of regular refinments. 

For \code{manifold@name=="sphere"}, we have the following grids:
there isn't exist \code{typegrid=="regular"} grid for a sphere, but
a \code{typegrid=="random"} grid (uniform density sample on the sphere of size \code{Ng*Ng}) and a \code{typegrid=="visualization"} grid (sphere-visualization grid on the sphere of size \code{6*Ng*Ng}, union
of the 6 domains centered around one of the 6 triply orthogonal poles, each domain are composed of the heights on the sphere (when they exists) corresponding to the regular mesh [-3/4,3/4]x[-3/4,3/4] of the others two cartesian coordinates).  

Finally, for \code{manifold@name=="hyperboloid"} we have:
no \code{typegrid=="regular"} grid on the hyperboloid, but a
\code{typegrid=="random"} grid (uniform density sample on the sphere of size \code{Ng*Ng}) and a 
\code{typegrid=="visualization"} grid (hyperboloid-vizualisation grid of size \code{Ng*Ng}, a domain of composed of the height
of the hyperboloid corresponding to the regular mesh [-3,3]x[-3,3] of the other two cartesian coordinates)
}

\author{
  Alexandre Brouste (\url{http://subaru.univ-lemans.fr/sciences/statist/pages_persos/Brouste/})
   and
  Sophie Lambert-Lacroix (\url{http://ljk.imag.fr/membres/Sophie.Lambert}).
}

\seealso{\code{\link{fieldsim}}.}

\examples{
# Load FieldSim library
library(FieldSim)

#Regular (or visualization) grid on plane
plane<-setManifold("plane")
plane@atlas<-constructgrid(plane,"regular",15)
str(plane@atlas)

#Random grid on the plane
plane@atlas<-constructgrid(plane,"random",15)
str(plane@atlas)

#For more examples of use of constructgrid, see fieldsim documentation.
}

