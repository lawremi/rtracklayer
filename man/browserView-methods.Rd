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
  browserView(object, range, track, ...)
}
\arguments{
  \item{object}{The object from which to get the views.}
  \item{range}{The \code{\link[GenomicRanges:GRanges-class]{GRanges}} or
    \code{\link[IRanges:RangesList-class]{RangesList}} to
    display. If there are multiple elements, a view is created for each
    element and a \code{\linkS4class{BrowserViewList}} is
    returned.} 
  \item{track}{List of track names to make visible in the view.}
  \item{\dots}{Arguments to pass to methods}
}
\section{Methods}{
  The following methods are defined by \pkg{rtracklayer}. 
  \describe{
    
    \item{object = "UCSCSession"}{ \code{browserView(object, range =
        range(object), track = trackNames(object), imagewidth = 800,
        ...)}: Creates a 
      \code{\linkS4class{BrowserView}} of \code{range} with visible
      tracks specified by \code{track}. The \code{imagewidth} parameter
      specifies the width of the track image in pixels. \code{track}
      may be an instance 
      of \code{\linkS4class{UCSCTrackModes}}. Arguments in \code{...}
      are passed to \code{\link{ucscTrackModes}} to create the
      \code{UCSCTrackModes} instance that will override modes
      indicated by the \code{track} parameter.
   }
    
}}
\examples{
\dontrun{
  session <- browserSession()
  browserView(session,
              GRangesForUCSCGenome("hg19", "chr2", IRanges(20000, 50000)))
  ## only view "knownGene" track
  browserView(session, track = "knownGene")
}
}
\keyword{methods}
