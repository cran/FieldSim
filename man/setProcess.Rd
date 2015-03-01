\encoding{latin1}

\name{setProcess}

\alias{setProcess}

\title{Construct usual processes on manifolds}

\usage{setProcess(name,parameter,values,manifold,covf)}

\arguments{

\item{name}{the name of the process (see details);}
\item{parameter}{the parameters of the process (see details);}
\item{values}{the values of the simulated (or given) sample path of the process;}
\item{manifold}{the manifold of which the process is defined;}
\item{covf}{the autocovariance function of the process.}
}

\value{an object of class \code{process}.}

\description{The function \code{setProcess} constructs usual processes on a specific manifold.}

\references{
A. Brouste, J. Istas and S. Lambert-Lacroix (2010) On simulation of manifold indexed fractional Gaussian fields. A. Brouste, J. Istas and S. Lambert-Lacroix (2014) Fractional Gaussian bridges with the package \code{FieldSim}.
}

\details{ 

We list here the different usual process. 

}

\author{Alexandre Brouste (\url{http://perso.univ-lemans.fr/~abrouste/})
   and
  Sophie Lambert-Lacroix (\url{http://membres-timc.imag.fr/Sophie.Lambert/}).
}

\seealso{\code{\link{fieldsim}}}

\examples{
# Load FieldSim library
library(FieldSim)

# Fractional Brownian field on [0,1]^2
plane.fBm<-setProcess("fBm-plane",0.7)
str(plane.fBm)

# Multifractional Brownian field on [0,1]^2
funcH<-function(xi){0.3+xi[1]*0.6}
plane.mBm<-setProcess("mBm-plane",funcH)
str(plane.mBm)

# Fractional Brownian sheet on [0,1]^2
#plane.fBs<-setProcess("fBs-plane",c(0.9,0.3))
#str(plane.fBs)

# Anisotropic fractional Brownian field on [0,1]^2
#plane.afBf<-setProcess("afBf-plane",list(H=0.7,theta1=pi/6,theta2=pi/3))
#str(plane.afBf)

# Bifractional fractional Brownian field on [0,1]^2
#plane.2pfBm<-setProcess("2pfBm-plane",list(H=0.7,K=0.5))
#str(plane.2pfBm)

# Spherical fractional Brownian field
#sphere.fBm<-setProcess("fBm-sphere",0.3)
#str(sphere.fBm)

# Fractional Brownian field on the hyperboloid
#hyper.fBm<-setProcess("fBm-hyperboloid",0.7)
#str(hyper.fBm)

# Bridge associated to the Fractional Brownian field on [0,1]^2
#Gamma<-matrix(c(1,0,0,0,1,1,1,1,1,1/2,1/2,0.5),3,4)
#bridge.plane.fBm<-setProcess("bridge-fBm-plane",list(Gamma=Gamma,par=0.9))
#str(bridge.plane.fBm)

# User defined process (see Brouste et al. 2010)

#sphere<-setManifold("sphere")
#user.sphere<-setProcess(name="user",manifold=sphere)

#parameter<-0.7
#acov<-function(xi,xj){exp(-#user.sphere@manifold@distance(xi,xj)^{2*user.sphere@parameter})}

#user.sphere@parameter<-parameter
#user.sphere@covf<-acov

#fieldsim(user.sphere)
#plot(user.sphere)
}

