\name{GraphTrackLine-class}
\docType{class}
\alias{GraphTrackLine-class}
\alias{coerce,character,GraphTrackLine-method}
\alias{coerce,GraphTrackLine,character-method}
\alias{coerce,GraphTrackLine,BasicTrackLine-method}
\alias{coerce,BasicTrackLine,GraphTrackLine-method}

\title{Class "GraphTrackLine"}
\description{ A \acronym{UCSC} track line for graphical tracks. }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("GraphTrackLine",
    ...)} or parsed from a character vector track line with
  \code{as(text, "GraphTrackLine")} or converted from a
  \code{\linkS4class{BasicTrackLine}} using \code{as(basic, "GraphTrackLine")}.
}
\section{Slots}{
  \describe{
    \item{\code{altColor}:}{Object of class \code{"integer"} giving an
      alternate color, as from \code{\link{col2rgb}}.}
    \item{\code{autoScale}:}{Object of class \code{"logical"} indicating
      whether to automatically scale to min/max of the data.}
    \item{\code{alwaysZero}:}{Object of class \code{"logical"} indicating
      whether to fix the lower limit of the Y axis at zero.}
    \item{\code{gridDefault}:}{Object of class \code{"logical"}
      indicating whether a grid should be drawn. }
    \item{\code{maxHeightPixels}:}{Object of class \code{"numeric"} of
      length three (max, default, min), giving the allowable range for
      the vertical height of the graph.}
    \item{\code{graphType}:}{Object of class \code{"character"},
      specifying the graph type, either "bar" or "points". }
    \item{\code{viewLimits}:}{Object of class \code{"numeric"} and of
      length two specifying the data range (min, max) shown in the graph.}
    \item{\code{yLineMark}:}{Object of class \code{"numeric"} giving the
      position of a horizontal line.}
    \item{\code{yLineOnOff}:}{Object of class \code{"logical"}
      indicating whether the \code{yLineMark} should be visible. }
    \item{\code{windowingFunction}:}{Object of class \code{"character"},
      one of "maximum", "mean", "minimum", for removing points when the
      graph shrinks.}
    \item{\code{smoothingWindow}:}{Object of class \code{"numeric"}
      giving the window size of a smoother to pass over the graph. }
    \item{\code{type}:}{Scalar \code{"character"} indicating the type of
      the track, either \code{"wig"} or \code{"bedGraph"}.}
    \item{\code{name}:}{Object of class \code{"character"} specifying
      the name of the track. }
    \item{\code{description}:}{Object of class \code{"character"}
      describing the track. }
    \item{\code{visibility}:}{Object of class \code{"character"}
      indicating the default visible mode of the track, see
      \code{\linkS4class{UCSCTrackModes}}.}
    \item{\code{color}:}{Object of class \code{"integer"} representing the
      track color (as from \code{\link{col2rgb}}).
    }
    \item{\code{priority}:}{Object of class \code{"numeric"} specifying
      the rank of this track. }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{TrackLine}"}, directly.
}
\section{Methods}{
  \describe{
    \item{as(object, "character")}{Export line to its string
      representation.}
    \item{as(object, \code{"\linkS4class{BasicTrackLine}"})}{Convert this
      line to a basic UCSC track line, using defaults for slots not held
      in common.}
  }
}
\references{ Official documentation:
  \url{http://genome.ucsc.edu/goldenPath/help/wiggle.html}. }
\author{ Michael Lawrence }
\seealso{
  \code{\link{export.wig}}, \code{\link{export.bedGraph}} for exporting
  graphical tracks.
}
\keyword{classes}
