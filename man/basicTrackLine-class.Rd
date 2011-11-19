\name{BasicTrackLine-class}
\docType{class}
\alias{BasicTrackLine-class}
\alias{coerce,character,BasicTrackLine-method}
\alias{coerce,BasicTrackLine,character-method}

\title{Class "BasicTrackLine"}
\description{The type of \acronym{UCSC} track line used to annotate most
  types of tracks (every type except Wiggle).}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("BasicTrackLine",
  ...)} or parsed from a character vector track line with \code{as(text,
  "BasicTrackLine")} or converted from a
\code{\linkS4class{GraphTrackLine}} using \code{as(wig, "BasicTrackLine")}.
}
\section{Slots}{
  \describe{
    \item{\code{itemRgb}:}{Object of class \code{"logical"} indicating
      whether each feature in a track uploaded as \acronym{BED} should
      be drawn in its specified color. }
    \item{\code{useScore}:}{Object of class \code{"logical"} indicating
      whether the data value should be mapped to color. }
    \item{\code{group}:}{Object of class \code{"character"} naming a
      group to which this track should belong. }
    \item{\code{db}:}{Object of class \code{"character"} indicating the
      associated genome assembly. }
    \item{\code{offset}:}{Object of class \code{"numeric"}, a number
      added to all positions in the track. }
    \item{\code{url}:}{Object of class \code{"character"} referring to
      additional information about this track. }
    \item{\code{htmlUrl}:}{Object of class \code{"character"} referring
      to an \acronym{HTML} page to be displayed with this track. }
    \item{\code{name}:}{Object of class \code{"character"} specifying
      the name of the track. }
    \item{\code{description}:}{Object of class \code{"character"}
      describing the track. }
    \item{\code{visibility}:}{Object of class \code{"character"}
      indicating the default visible mode of the track, see
      \code{\linkS4class{UCSCTrackModes}}.}
    \item{\code{color}:}{Object of class \code{"integer"} representing the
      track color (as from \code{\link{col2rgb}}). }
    \item{\code{colorByStrand}:}{Object of class \code{"matrix"}
       with two columns, as from \code{col2rgb}. The two colors
       indicate the color for each strand (positive, negative). }
    \item{\code{priority}:}{Object of class \code{"numeric"} specifying
      the rank of the track. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{TrackLine}"}, directly.
}
\section{Methods}{
  \describe{
    \item{as(object, "character")}{Export line to its string
      representation.}
    \item{as(object, \code{"\linkS4class{GraphTrackLine}"})}{Convert this
      line to a graph track line, using defaults for slots not held
      in common.}
  }
}
\references{
  \url{http://genome.ucsc.edu/goldenPath/help/customTrack.html#TRACK}
  for the official documentation. }
\author{ Michael Lawrence }

\seealso{
  \code{\linkS4class{GraphTrackLine}} for Wiggle/bedGraph tracks.
}
\keyword{classes}
