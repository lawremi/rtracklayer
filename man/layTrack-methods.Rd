\name{layTrack-methods}
\docType{methods}
\alias{layTrack}
\alias{layTrack-methods}
\alias{layTrack,BrowserSession,RangedData-method}
\alias{layTrack,BrowserSession,RangedDataList-method}
\alias{layTrack,UCSCSession,RangedDataList-method}
\title{ Laying tracks }
\description{
  Methods for loading \code{\linkS4class{RangedData}} instances (tracks)
  into genome browsers.
}
\usage{
layTrack(object, track, name = deparse(substitute(track)), view = FALSE, ...)
}
\arguments{
  \item{object}{A \code{\linkS4class{BrowserSession}} into which the
    track is loaded.}
  \item{track}{The track(s) to load.}
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

    \item{object = "BrowserSession", track = "RangedData"}{ Load this
      track into the session. }

    \item{object = "BrowserSession", track = "RangedDataList"}{ Load all tracks
      into the session. }

    \item{object = "UCSCSession", track = "RangedDataList"}{
      \code{layTrack(object, track, name = deparse(substitute(track)),
        view = FALSE, format = "gff")}:
      Load the tracks into the session using the specified
      \code{format}.
    }
  }
}
\seealso{\code{\link{track}} for getting a track from a session.}
\examples{
\dontrun{
  session <- browserSession()
  track <- import(system.file("tests", "v1.gff", package = "rtracklayer"))
  layTrack(session, track, "My Track")
}
}
\keyword{methods}
