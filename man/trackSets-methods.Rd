\name{trackSets-methods}
\docType{methods}
\alias{trackSets}
\alias{trackSets-methods}
\alias{trackSets,trackSet-method}
\alias{trackSets,missing-method}
\alias{trackSets,list-method}
\title{Getting a list of trackSets}
\description{
 Methods for creating and getting a \code{\linkS4class{trackSets}}
 object: a list of \code{\linkS4class{trackSet}} instances.
}
\section{Methods}{
  The following methods are defined by \pkg{rtracklayer} for the
  \code{trackSets(object, ...)} generic.
  
  \describe{
    \item{object = "trackSet"}{ Create an instance with track
      \code{object} and any additional tracks in \code{\dots}. }
    \item{object = "missing"}{ Creates an empty instance. }
    \item{object = "list"}{ Coerces a regular list of \code{trackSet}
      instances. }
  }
}
\keyword{methods}
