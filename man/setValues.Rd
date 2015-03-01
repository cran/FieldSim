\encoding{latin1}

\name{setValues}
\alias{setValues}
\alias{setValues,process-method}

\title{Set the values of an object of class \code{process}}
\usage{setValues(process,values)}

\arguments{
\item{process}{an S4 object process;}
\item{values}{the values of the process on the atlas.}
}

\description{The function \code{setValues} set the values of the process.}

\references{A. Brouste, J. Istas and S. Lambert-Lacroix (2015). Fractional Gaussian bridges with the package \code{FieldSim}.}

\details{Statistical tools developed in the Fieldsim package allows real dataset inserted in the model process. Consequently the user can set the values of the process. Parameter will be forgotten using statistical command. This values will be erase with the use of fieldsim.
}

\author{Alexandre Brouste (\url{http://perso.univ-lemans.fr/~abrouste/}) and
 Sophie Lambert-Lacroix (\url{http://membres-timc.imag.fr/Sophie.Lambert/}).}

\seealso{\code{\link{fieldsim}}.}

\examples{
# Load FieldSim library
library(FieldSim)

#Dataset (to do)
plane.fBm<-setProcess("fBm-plane",0.6)
fieldsim(plane.fBm)
sample<-plane.fBm@values

plane.fBm.2<-setProcess("fBm-plane",0.7)
setValues(plane.fBm.2,sample)
}

