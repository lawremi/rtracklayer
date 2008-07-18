\name{ucscTrackSet-class}
\docType{class}
\alias{ucscTrackSet-class}
\alias{export.bed,ucscTrackSet-method}
\alias{export.gff,ucscTrackSet-method}
\alias{export.ucsc,ucscTrackSet-method}

\title{Class "ucscTrackSet"}
\description{Each track in \acronym{UCSC} has an associated
  \code{\linkS4class{ucscTrackLine}} that contains metadata on the track.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ucscTrackSet",
    assayData, phenoData, featureData, experimentData, annotation,
    ...)}. See \code{\linkS4class{trackSet}} for more details.
}
\section{Slots}{
  \describe{
    \item{\code{trackLine}:}{Object of class \code{"ucscTrackLine"}
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
\section{Extends}{
Class \code{"\linkS4class{trackSet}"}, directly.
Class \code{"\link[Biobase:eSet-class]{eSet}"}, by class "trackSet", distance 2.
Class \code{"\link[Biobase:VersionedBiobase-class]{VersionedBiobase}"}, by class "trackSet", distance 3.
Class \code{"\link[Biobase:Versioned-class]{Versioned}"}, by class "trackSet", distance 4.
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
