\name{tracks-methods}
\docType{methods}
\alias{trackNames}
\alias{trackNames<-}
\alias{trackNames<--methods}
\alias{trackNames<-,UCSCTrackModes-method}
\alias{trackNames<-,UCSCView-method}
\alias{trackNames-methods}
\alias{trackNames,UCSCSession-method}
\alias{trackNames,UCSCTrackModes-method}
\alias{trackNames,UCSCView-method}
\title{ Accessing track names }
\description{
  Methods for getting and setting track names.
}
\section{Methods}{
\describe{
  The following methods are defined by \pkg{rtracklayer} for \bold{getting}
  track names via the generic \code{trackNames(object, ...)}.
  
  \item{object = "UCSCSession"}{ Get the tracks loaded in the session. }

  \item{object = "UCSCTrackModes"}{ Get the visible tracks according to
    the modes (all tracks not set to "hide").}

  \item{object = "UCSCView"}{ Get the visible tracks in the view. }

  The following methods are defined by \pkg{rtracklayer} for \bold{setting}
  track names via the generic \code{trackNames(object) <- value}.
  
  \item{object = "UCSCTrackModes"}{ Sets the tracks that should be
    visible in the modes. All specified tracks with mode "hide" in
    \code{object} are set to mode "full". Any tracks in \code{object} that
    are not specified in the value are set to "hide". No other modes are
    changed.}

  \item{object = "UCSCView"}{ Sets the visible tracks in the view. This
    opens a new web browser with only the specified tracks visible. }
}}
\keyword{methods}