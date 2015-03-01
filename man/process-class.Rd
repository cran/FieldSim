\name{process-class}
\docType{class}

\alias{process}
\alias{process-class}

\title{process class}
\description{
  The \code{process} class is a class of the  \pkg{FieldSim} package.  
}
\section{Slots}{
  \describe{
 
\item{\code{name}:}{is the name of the manifold (a character string).}

\item{\code{values}:}{the values of the simulated (or given) sample path of the process ().} 

\item{\code{parameter}:}{is the origin fixed on the manifold (a matrix)} 

\item{\code{manifold}:}{is the distance set on the manifold (a function).} 
 
\item{\code{covf}:}{is the mesh (a matrix).}


}	 
}

\details{

Several names for slot \code{name} are reserved for classical fractional Gaussian processes: \code{"fBm"} for fractional Brownian motion,  \code{"mBm"} for multifractional Brownian motion, \code{"2pfBm"} for the standard bi-fractional Brownian motion, \code{"stdfBm"} for the space-time deformed fractional Brownian motion, \code{"afBf"} for anisotropic fractional Brownian field, \code{"fBs"} for fractional Brownian sheet  and \code{"bridge"} for all kind of bridges.  

The slot \code{manifold} contains an object of class \code{manifold} (see \link{manifold-class}).

The slot \code{parameter} that contains all the parameter associated to the covariance function of the process. Here are the classical parameter associated to classical processes. For instance, \code{"fbm"} has parameter \code{numeric}, \code{"mbm"} has parameter \code{function}, \code{"2pfbm"} has parameter \code{list(H=numeric, K=numeric)}, 
\code{"stdfbm"} has parameter \code{list(H=numeric, sigma=function, tau=function)}, 
\code{"afBf"} has parameter \code{list(H=numeric, theta1=numeric, theta2=numeric)},
\code{"fBs"} has parameter \code{vector} and 
\code{"bridge"} has  \code{list(Gamma=matrix, R=function, Tp=..., par=list(...))}.
}

\seealso{\code{\link{setProcess}},\code{\link{manifold-class}}.}

\author{Alexandre Brouste}
\keyword{classes}
