\name{genomeSequence-methods}
\docType{methods}
\alias{genomeSequence}
\alias{genomeSequence-methods}
\alias{genomeSequence,UCSCSession-method}
\title{ Retrieving a genome sequence }
\description{
  Methods for retrieving the sequence of a
  \code{\linkS4class{RangesList}} from an object.
}
\section{Methods}{
  The following methods are defined by \pkg{rtracklayer} for
  \code{genomeSequence(object, range = range(object), ...)}.
  \describe{
    \item{object = "UCSCSession"}{
      \code{genomeSequence(object,
        range = range(object), track = "Assembly")}:
      Gets the sequence in \code{range} and \code{track} from the
      session. 
    }
}}
\seealso{\code{\link{laySequence}} for storing sequences.}
\keyword{methods}
