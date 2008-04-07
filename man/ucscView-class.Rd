\name{ucscView-class}
\docType{class}
\alias{ucscView-class}

\title{Class "ucscView"}
\description{An object representing a view of a genome in the
  \acronym{UCSC} browser.}
\section{Objects from the Class}{
  Objects are created by calling \code{\link{browserView}(session, ...)}
  where \code{session} is a \code{\linkS4class{ucscSession}}. 
}
\section{Slots}{
  \describe{
    \item{\code{hgsid}:}{Object of class \code{"numeric"}, which
      identifies this view to \acronym{UCSC}.}
    \item{\code{session}:}{Object of class \code{"browserSession"} to
      which this view belongs.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{browserView}"}, directly.
}
\section{Methods}{
  \describe{
    \item{\code{\link{activeView}(object)}}{ Obtains a logical
      indicating whether this view is the active view. }
    \item{\code{\link{genomeSegment}(object)}}{ Obtains the
      \code{\linkS4class{genomeSegment}} displayed by this view. }
    \item{\code{\link{genomeSegment}(object) <- value}}{ Sets the
      \code{\linkS4class{genomeSegment}} displayed by this view. }
    \item{\code{\link{tracks}(object)}}{ Gets the names of the visible
      tracks in this view. }
    \item{\code{\link{tracks}(object) <- value}}{ Sets the visible
      tracks by name. }
    \item{\code{\link{ucscTrackModes}(object)}}{ Obtains the
      \code{\linkS4class{ucscTrackModes}} for this view. }
    \item{\code{\link{ucscTrackModes}(object) <- value}}{ Sets the
      \code{\linkS4class{ucscTrackModes}} for this view. The
      \code{value} may be either a \code{ucscTrackModes} instance or a
      character vector that will be coerced by a call to
      \code{\link{ucscTrackModes}}. }
  }
}
\author{ Michael Lawrence }
\seealso{
  \code{\link{browserView}} for creating instances of this class.
}
\keyword{classes}
