\name{genomeSequence-methods}
\docType{methods}
\alias{genomeSequence}
\alias{genomeSequence-methods}
\alias{genomeSequence,ucscSession-method}
\title{ Retrieving a genome sequence }
\description{
  Methods for retrieving the sequence of a
  \code{\linkS4class{genomeSegment}} from an object.
}
\section{Methods}{
  The following methods are defined by \pkg{rtracklayer} for
  \code{genomeSequence(object, segment = genomeSegment(object), ...)}.
  \describe{
    \item{object = "ucscSession"}{
      \code{genomeSequence(object,
        segment = genomeSegment(object), track = "Assembly")}:
      Gets the sequence in \code{segment} and \code{track} from the
      session. 
    }
}}
\seealso{\code{\link{laySequence}} for storing sequences.}
\keyword{methods}
