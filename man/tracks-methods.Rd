\name{tracks-methods}
\docType{methods}
\alias{tracks}
\alias{tracks<-}
\alias{tracks<--methods}
\alias{tracks<-,ucscTrackModes-method}
\alias{tracks<-,ucscView-method}
\alias{tracks-methods}
\alias{tracks,ucscSession-method}
\alias{tracks,ucscTrackModes-method}
\alias{tracks,ucscView-method}
\title{ Accessing track names }
\description{
  Methods for getting and setting track names.
}
\section{Methods}{
\describe{
  The following methods are defined by \pkg{rtracklayer} for \bold{getting}
  track names via the generic \code{tracks(object, ...)}.
  
  \item{object = "ucscSession"}{ Get the tracks loaded in the session. }

  \item{object = "ucscTrackModes"}{ Get the visible tracks according to
    the modes (all tracks not set to "hide").}

  \item{object = "ucscView"}{ Get the visible tracks in the view. }

  The following methods are defined by \pkg{rtracklayer} for \bold{setting}
  track names via the generic \code{tracks(object) <- value}.
  
  \item{object = "ucscTrackModes"}{ Sets the tracks that should be
    visible in the modes. All specified tracks with mode "hide" in
    \code{object} are set to mode "full". Any tracks in \code{object} that
    are not specified in the value are set to "hide". No other modes are
    changed.}

  \item{object = "ucscView"}{ Sets the visible tracks in the view. This
    opens a new web browser with only the specified tracks visible. }
}}
\keyword{methods}