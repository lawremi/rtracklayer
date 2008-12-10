\name{UCSCData-class}
\docType{class}
\alias{UCSCData-class}
\alias{export.bed,UCSCData-method}
\alias{export.gff,UCSCData-method}
\alias{export.ucsc,UCSCData-method}
\alias{initialize,UCSCData-method}
\alias{show,UCSCData-method}

\title{Class "UCSCData"}
\description{Each track in \acronym{UCSC} has an associated
  \code{\linkS4class{TrackLine}} that contains metadata on the track.}
\section{Constructor}{
}
\section{Slots}{
  \describe{
    \item{\code{trackLine}:}{Object of class \code{"TrackLine"}
      holding track metadata. }
    \item{\code{genome}:}{Object of class \code{"character"} identifying
      the genome of this track, e.g. "hg18".}
    \item{\code{assayData}:}{Object of class \code{"AssayData"} holding
      the experimental measurements under \code{dataVals}. }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}
      holding the experimental design matrix. }
    \item{\code{featureData}:}{Object of class
      \code{"AnnotatedDataFrame"} holding feature information, generally
      including columns \code{chrom} (chromosome), \code{start}
      (start position), \code{end} (end position) and
      \code{strand} (strand on the DNA: "+", "-", or \code{NA}).}
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} holding
    experimental metadata.} 
    \item{\code{annotation}:}{Object of class \code{"character"}
      referring to the annotation dataset. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}
      holding version information. }
  }
}
\section{Methods}{
  \describe{
    \item{\code{\link{export.bed}(object, con, wig = FALSE, trackLine =
        !wig)}}{ Exports the track and its track line (if
      \code{trackLine} is \code{TRUE}) to \code{con} in the Browser
      Extended Display (\acronym{BED}) format. }
    \item{\code{\link{export.gff}(object)}}{ Exports
    the track and its track line (as a comment) to \code{con} in the
    General Feature Format (\acronym{GFF}). }
    \item{\code{\link{export.ucsc}(object)}}{
      Exports the track and its track line to \code{con} in the
      \acronym{UCSC} meta-format.  }
  }
}
\author{ Michael Lawrence }
\seealso{
  \code{\link{import}} and \code{\link{export}} for reading and writing
  tracks to and from connections (files), respectively.
}
\keyword{classes}
