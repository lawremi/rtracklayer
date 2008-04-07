\name{trackSets-class}
\docType{class}
\alias{trackSets-class}

\title{Class "trackSets"}
\description{A list of \code{\linkS4class{trackSet}} instances.}
\section{Objects from the Class}{
  Objects are normally created using the \code{\link{trackSets}} function.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} containing the tracks. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\link[Biobase:AssayData-class]{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{\code{\link{export.ucsc}(object, con,
        subformat = c("gff", "wig", "bed"), name = names(object), ...)}}{Export
      all tracks in \acronym{UCSC} format. }
    \item{\code{\link{genomeSegment}(object)}}{Obtain
      the union of the segments spanned by each track in the list. }
  }
}
\author{ Michael Lawrence }

\seealso{
  \code{\link{trackSets}} for obtaining instances of this class.

  \code{\link{layTrack}} for loading the tracks into a genome browser.
}
\keyword{classes}
