\name{track<--methods}
\docType{methods}
\alias{track<-}
\alias{track<--methods}
\alias{track<-,BrowserSession,ANY-method}
\title{ Laying tracks }
\description{
  Methods for loading tracks into genome browsers.
}
\usage{
track(object, ...) <- value
}
\arguments{
  \item{object}{A \code{\linkS4class{BrowserSession}} into which the
    track is loaded.}
  \item{value}{The track(s) to load.}
  \item{\dots}{Arguments to pass on to methods. Can be:
    \itemize{
      \item{name}{The name(s) of the track(s) being loaded.}
      \item{view}{Whether to create a view of the track after loading it.}
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
