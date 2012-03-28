\name{UCSCData-class}
\docType{class}
\alias{UCSCData-class}
\alias{export.bed,UCSCData,characterORconnection-method}
\alias{export.bed15,UCSCData-method}
\alias{export.gff,UCSCData,characterORconnection-method}
\alias{export.ucsc,UCSCData,characterORconnection-method}
\alias{initialize,UCSCData-method}
\alias{show,UCSCData-method}
\alias{coerce,RangedData,UCSCData-method}
\alias{coerce,UCSCData,GRanges-method}

\title{Class "UCSCData"}
\description{Each track in \acronym{UCSC} has an associated
  \code{\linkS4class{TrackLine}} that contains metadata on the track.}
\section{Slots}{
  \describe{
    \item{\code{trackLine}:}{Object of class \code{"TrackLine"}
      holding track metadata. }
  }
}
\section{Methods}{
  \describe{
    \item{\code{\link{export.bed}(object, con, variant = c("base",
        "bedGraph", "bed15"), color, trackLine = TRUE, ...)}}{
      Exports the track and its track line (if
      \code{trackLine} is \code{TRUE}) to \code{con} in the Browser
      Extended Display (\acronym{BED}) format. The arguments
      in \code{...} are passed to \code{\link{export.ucsc}}.
    }
    \item{\code{\link{export.bed15}(object, con, expNames = NULL, ...)}}{
      Exports the track and its track line (if
      \code{trackLine} is \code{TRUE}) to \code{con} in the Bed15
      format. The data is taken from the columns named in
      \code{expNames}, which defaults to the \code{expNames} in the
      track line, if any, otherwise all column names. The arguments
      in \code{...} are passed to \code{\link{export.ucsc}}.
    }
    \item{\code{\link{export.gff}(object)}}{
      Exports the track and its track line (as a comment) to \code{con} in the
      General Feature Format (\acronym{GFF}).
    }
    \item{\code{\link{export.ucsc}(object, con, subformat, ...)}}{
      Exports the track and its track line to \code{con} in the
      \acronym{UCSC} meta-format.  }
    \item{\code{as(object, "UCSCData")}}{Constructs a \code{UCSCData}
      from a \code{RangedData} instance, by adding a default track line
      and ensuring that the sequence/chromosome names are compliant with
      UCSC conventions. If there is a numeric score, the track line type
      is either "bedGraph" or "wig", depending on the feature
      density. Otherwise, "bed" is chosen.
    }
  }
}
\author{ Michael Lawrence }
\seealso{
  \code{\link{import}} and \code{\link{export}} for reading and writing
  tracks to and from connections (files), respectively.
}
\keyword{classes}
