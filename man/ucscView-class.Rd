\name{UCSCView-class}
\docType{class}
\alias{UCSCView-class}
\alias{range,UCSCView-method}
\alias{range<-,UCSCView-method}
\alias{visible,UCSCView-method}
\alias{visible<-,UCSCView-method}

\title{Class "UCSCView"}
\description{An object representing a view of a genome in the
  \acronym{UCSC} browser.}
\section{Objects from the Class}{
  Calling \code{\link{browserView}(session, range = range(object), track
    = trackNames(object), ...)}
  creates \code{\linkS4class{BrowserView}} of \code{range} with visible
  tracks specified by \code{track}. \code{track} may be an instance of
  \code{\linkS4class{UCSCTrackModes}}. Arguments in \code{...} should
  match parameters to a
  \code{\link{ucscTrackModes}} method for creating a
  \code{UCSCTrackModes} instance that will be merged with and override
    modes indicated by the \code{track} parameter.
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
      \code{\link[GenomicRanges:GRanges-class]{GRanges}} displayed by
      this view. } 
    \item{\code{\link{range}(object) <- value}}{ Sets the
      \code{\link[GenomicRanges:GRanges-class]{GRanges}} or
      \code{\link[IRanges:RangesList-class]{RangesList}} displayed by
      this view.  } 
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
