\name{dataVals-methods}
\docType{methods}
\alias{dataVals}
\alias{dataVals-methods}
\alias{dataVals,trackSet-method}
\title{ Get data values }
\description{
  Methods for getting data values.
}
\section{Methods}{
  The following methods are defined by \pkg{rtracklayer}. 
\describe{
  \item{object = "trackSet"}{ Get the experimental measurements from a track.}
}}
\examples{
  track <- import(system.file("tests", "bed.wig", package = "rtracklayer"))
  dataVals(track)
}
\keyword{methods}