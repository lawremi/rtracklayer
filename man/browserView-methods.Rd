\name{browserView-methods}
\docType{methods}
\alias{browserView}
\alias{browserView-methods}
\alias{browserView,UCSCSession-method}
\title{Getting browser views}
\description{
  Methods for creating and getting browser views.
}
\usage{
  browserView(object, range = range(object),
              track = trackNames(object), ...)
}
\arguments{
  \item{object}{The object from which to get the views.}
  \item{range}{The \code{\linkS4class{RangesList}} to display.}
  \item{track}{List of track names to make visible in the view.}
  \item{\dots}{Arguments to pass to methods}
}
\section{Methods}{
  The following methods are defined by \pkg{rtracklayer}. 
  \describe{
    
    \item{object = "UCSCSession"}{ \code{browserView(object, range =
        range(object), track = trackNames(object), ...)}: Creates a
      \code{\linkS4class{BrowserView}} of \code{range} with visible
      tracks specified by \code{track}. \code{track} may be an instance
      of \code{\linkS4class{UCSCTrackModes}}. Arguments in \code{...} should
      override slots in \code{range} or else match parameters
      to a \code{\link{ucscTrackModes}} method for creating a
      \code{UCSCTrackModes} instance that will override modes indicated
      by the \code{track} parameter. }
    
}}
\examples{
\dontrun{
  session <- browserSession()
  browserView(session, GenomicRanges(start = 20000, end = 50000,
                                     range = range(session)))
  ## equivalent to above, but shorter
  browserView(session, start = 20000, end = 50000)
  ## only view "knownGene" track
  browserView(session, track = "knownGene")
}
}
\keyword{methods}
