\name{process-method}
\docType{class}

\alias{initialize, process-method}

\title{process method}
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


\seealso{\code{\link{setProcess}},\code{\link{manifold-class}}.}

\author{Alexandre Brouste}
\keyword{classes}
