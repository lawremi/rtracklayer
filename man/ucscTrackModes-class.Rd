\name{UCSCTrackModes-class}
\docType{class}
\alias{UCSCTrackModes-class}
\alias{[,UCSCTrackModes,ANY,ANY-method}
\alias{[<-,UCSCTrackModes,ANY,ANY,ANY-method}

\title{Class "UCSCTrackModes"}
\description{A vector of view modes ("hide", "dense", "full",
  "pack", "squish") for each track in a \acronym{UCSC} view.}
\section{Objects from the Class}{
  Objects may be created by calls of the form
  \code{\link{ucscTrackModes}(object =
    character(), hide = character(), dense = character(), pack = character(),
    squish = character(), full = character())},
  where \code{object} should be a character vector of mode names (with
  its \code{names} attribute specifying the corresponding track names).
  The other parameters should contain track names that override the
  modes in \code{object}. Later parameters override earlier ones, so,
  for example, if a track is named in \code{hide} and \code{full},
  it is shown in the full view mode.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"character"} holding the
      modes ("hide", "dense", "full", "pack", "squish"), with its
      \code{names} attribute holding corresponding track names.}
    \item{\code{labels}:}{Object of class \code{"character"} holding
      labels (human-readable names) corresponding to each track/mode.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{character}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "character", distance 2.
}

\section{Methods}{
  \describe{
    \item{\code{\link{trackNames}(object)}}{ Gets the names of the visible
      tracks (those that do not have mode "hide"). }
    \item{\code{\link{trackNames}(object) <- value}}{ Sets the names of the
      visible tracks. Any tracks named in \code{value} are set to "full"
      if the are currently set to "hide" in this object. Any tracks not
      in \code{value} are set to "hide". All other modes are preserved.}
    \item{\code{object[i]}}{Gets the track mode of the tracks indexed
      by \code{i}, which can be any type of index supported by character
      vector subsetting. If \code{i} is a character vector, it indexes
      first by the internal track IDs (the \code{names} on \code{.Data})
      and then by the user-level track names (the \code{labels} slot).}
    \item{\code{object[i] <- value}}{Sets the track modes indexed by
      \code{i} (in the same way as in \code{object[i]} above) to those
      specified in \code{value}.} 
  }
}
\author{ Michael Lawrence }
\seealso{
  \code{\linkS4class{UCSCView}} on which track view modes may be set.
}
\keyword{classes}
