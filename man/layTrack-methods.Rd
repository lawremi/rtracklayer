\name{track<--methods}
\docType{methods}
\alias{track<-}
\alias{track<--methods}
\alias{track<-,BrowserSession,RangedData-method}
\alias{track<-,BrowserSession,ANY-method}
\alias{track<-,BrowserSession,RangedDataList-method}
\alias{track<-,UCSCSession,RangedDataList-method}
\title{ Laying tracks }
\description{
  Methods for loading \code{\link[IRanges:RangedData-class]{RangedData}} instances (tracks)
  into genome browsers.
}
\usage{
\S4method{track}{BrowserSession,RangedData}(object, name = deparse(substitute(track)), view = FALSE, ...) <- value
}
\arguments{
  \item{object}{A \code{\linkS4class{BrowserSession}} into which the
    track is loaded.}
  \item{value}{The track(s) to load.}
  \item{name}{The name(s) of the track(s) being loaded.}
  \item{view}{Whether to create a view of the track after loading it.}
  \item{\dots}{Arguments to pass on to methods.}
}
\section{Methods}{
  The following methods are defined by \pkg{rtracklayer}. A browser
  session implementation must implement a method for either
  \code{RangedData} or \code{RangedDataList}. The base \code{browserSession}
  class will delegate appropriately.
  
  \describe{

    \item{object = "BrowserSession", value = "RangedData"}{ Load this
      track into the session. }

    \item{object = "BrowserSession", value = "RangedDataList"}{ Load all tracks
      into the session. }

    \item{object = "UCSCSession", value = "RangedDataList"}{
      \code{track(object, name = deparse(substitute(track)),
        view = FALSE, format = "gff", ...) <- value}:
      Load the tracks into the session using the specified
      \code{format}. The arguments in \code{...} are passed on to
      \code{\link{export.ucsc}}, so they could be slots in a
      \code{\linkS4class{TrackLine}} subclass or parameters to pass
      on to the export function for \code{format}.
    }
  }
}
\seealso{\code{\link{track}} for getting a track from a session.}
\examples{
\dontrun{
  session <- browserSession()
  track <- import(system.file("tests", "v1.gff", package = "rtracklayer"))
  track(session, "My Track") <- track
}
}
\keyword{methods}
