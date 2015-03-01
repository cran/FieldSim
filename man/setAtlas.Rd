\encoding{latin1}
\name{setAtlas}
\alias{setAtlas}
\alias{setAtlas,manifold-method}
\alias{setAtlas,process-method}
\title{Construct usual grids on some specific manifolds}
\usage{setAtlas(object,gridtype,Ng)}
\arguments{
\item{object}{a object of class process or manifold;}
\item{gridtype}{the type of the grid, possible choice \code{"regular"}, \code{"random"} or \code{"visualization"};}
\item{Ng}{parameter of the size of the grid, see details.}}

\description{The function \code{setAtlas} constructs usual grids on manifold.}

\references{
A. Brouste, J. Istas and S. Lambert-Lacroix (2010). On simulation of manifold indexed fractional Gaussian fields.}

\details{ 
We list here the different implemented grids. For 
\code{manifold@name=="plane"} we have the \code{gridtype=="regular"} grid (with the parameter \code{Ng} returns a regular grid on [0,1]x[0,1] of size \code{Ng*Ng}), the \code{gridtype=="random"} grid (uniform random choice of the both coordinates on [0,1], grid of size \code{Ng*Ng}) and the \code{gridtype=="visualization"} grid, of size \code{ (2^Ng +1)*(2^Ng +1)} composed of regular refinments. 

For \code{manifold@name=="sphere"}, we have the following grids:
there isn't exist \code{gridtype=="regular"} grid for a sphere, but
a \code{gridtype=="random"} grid (uniform density sample on the sphere of size \code{Ng*Ng}) and a \code{gridtype=="visualization"} grid (sphere-visualization grid on the sphere of size \code{6*Ng*Ng}, union
of the 6 domains centered around one of the 6 triply orthogonal poles, each domain are composed of the heights on the sphere (when they exists) corresponding to the regular mesh [-3/4,3/4]x[-3/4,3/4] of the others two cartesian coordinates).  

Finally, for \code{manifold@name=="hyperboloid"} we have:
no \code{gridtype=="regular"} grid on the hyperboloid, but a
\code{gridtype=="random"} grid (uniform density sample on the sphere of size \code{Ng*Ng}) and a 
\code{gridtype=="visualization"} grid (hyperboloid-vizualisation grid of size \code{Ng*Ng}, a domain of composed of the height
of the hyperboloid corresponding to the regular mesh [-3,3]x[-3,3] of the other two cartesian coordinates)
}

\author{Alexandre Brouste (\url{http://perso.univ-lemans.fr/~abrouste/})
   and
  Sophie Lambert-Lacroix (\url{http://membres-timc.imag.fr/Sophie.Lambert/}).
}

\seealso{\code{\link{fieldsim}}.}

\examples{
# Load FieldSim library
library(FieldSim)

# Example on the line manifold

line<-setManifold("line")
str(line)

setAtlas(line,"regular",200)
str(line)

#Example on the fractional Brownian motion

line.fBm<-setProcess("fBm-line",0.7)
str(line.fBm)

setAtlas(line.fBm,"regular",200)
str(line.fBm)

setAtlas(line.fBm,"random",100)
str(line.fBm)

setAtlas(line.fBm,"finer",9)
str(line.fBm)

setAtlas(line.fBm,"visualization",9)
str(line.fBm)

}

