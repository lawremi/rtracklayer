\name{chrid-methods}
\docType{methods}
\alias{chrid}
\alias{chrid-methods}
\alias{chrid,trackSet-method}
\alias{chrid,character-method}
\alias{chrid,ANY-method}
\title{ Getting chromosome ids }
\description{
 Methods for creating and getting chromosome ids.
}
\section{Methods}{
  The following methods are defined by \pkg{rtracklayer}. 
  \describe{
  
    \item{object = "trackSet"}{ Gets the chromosomes of the features as a
      \code{\linkS4class{chrid}} instance. }
    
    \item{object = "character"}{ Coerces vector of string chromosome
      identifiers to a \code{\linkS4class{chrid}} instance. }

    \item{object = "ANY"}{ Attempts to coerce \code{object} by first
      coercing it to a character vector. }
}}
\seealso{Subset \code{\linkS4class{trackSet}} instances by chromosome
  identifier with \code{\link{[.trackSet}}}

\examples{
  chrid("chr22")
  track <- import(system.file("tests", "bed.wig", package = "rtracklayer"))
  chrid(track)
}
\keyword{methods}

