\name{UCSCView-class}
\docType{class}
\alias{UCSCView-class}
\alias{range,UCSCView-method}
\alias{range<-,UCSCView,RangesList-method}
\alias{visible,UCSCView-method}
\alias{visible<-,UCSCView-method}

\title{Class "UCSCView"}
\description{An object representing a view of a genome in the
  \acronym{UCSC} browser.}
\section{Objects from the Class}{
  Objects are created by calling \code{\link{browserView}(session, ...)}
  where \code{session} is a \code{\linkS4class{UCSCSession}}. 
}
\section{Slots}{
  \describe{
    \item{\code{hgsid}:}{Object of class \code{"numeric"}, which
      identifies this view to \acronym{UCSC}.}
    \item{\code{session}:}{Object of class \code{"BrowserSession"} to
      which this view belongs.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{BrowserView}"}, directly.
}
\section{Methods}{
  \describe{
    \item{\code{\link{activeView}(object)}}{ Obtains a logical
      indicating whether this view is the active view. }
    \item{\code{\link{range}(object)}}{ Obtains the
      \code{\link[IRanges:RangesList-class]{RangesList}} displayed by
      this view. } 
    \item{\code{\link{range}(object) <- value}}{ Sets the
      \code{\link[IRanges:RangesList-class]{RangesList}} displayed by
      this view. } 
    \item{\code{\link{trackNames}(object)}}{ Gets the names of the visible
      tracks in this view. }
    \item{\code{\link{trackNames}(object) <- value}}{ Sets the visible
      tracks by name. }
    \item{\code{visible(object)}}{Get a named logical vector indicating
      whether each track is visible.
    }
    \item{\code{visible(object) <- value}}{Set a logical vector
      indicating the visibility of each track, in the same order as
      returned by \code{visible(object)}.
    }
    \item{\code{\link{ucscTrackModes}(object)}}{ Obtains the
      \code{\linkS4class{UCSCTrackModes}} for this view. }
    \item{\code{\link{ucscTrackModes}(object) <- value}}{ Sets the
      \code{\linkS4class{UCSCTrackModes}} for this view. The
      \code{value} may be either a \code{UCSCTrackModes} instance or a
      character vector that will be coerced by a call to
      \code{\link{ucscTrackModes}}. }
    
  }
}
\author{ Michael Lawrence }
\seealso{
  \code{\link{browserView}} for creating instances of this class.
}
\keyword{classes}
