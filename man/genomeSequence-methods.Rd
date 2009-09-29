\name{getSeq-methods}
\docType{methods}
\alias{getSeq,UCSCSession-method}
\title{ Retrieving a genome sequence }
\description{
  Methods for retrieving the sequence of a
  \code{\link[IRanges:RangesList-class]{RangesList}} from an object.
}
\section{Methods}{
  The following methods are defined by \pkg{rtracklayer} for
  \code{getSeq(object, range = range(object), ...)}.
  \describe{
    \item{object = "UCSCSession"}{
      \code{getSeq(object,
        range = range(object), track = "Assembly")}:
      Gets the sequence in \code{range} and \code{track} from the
      session. 
    }
}}
\seealso{\code{\link{sequence<-}} for storing sequences.}
\keyword{methods}
