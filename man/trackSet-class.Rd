\name{trackSet-class}
\docType{class}
\alias{trackSet-class}
\alias{[.trackSet}
\alias{[,trackSet-method}
\alias{coerce,trackSet,data.frame-method}
\alias{initialize,trackSet-method}

\alias{chrom}
\alias{strand}
\alias{chrom,trackSet-method}
\alias{start,trackSet-method}
\alias{end,trackSet-method}
\alias{strand,trackSet-method}

\title{Class "trackSet"}
\description{ An object representing a genome annotation track. Based on
  \code{\link[Biobase:eSet-class]{eSet}}, with feature information
  (chrom, start, end, strand, ...) stored in the \code{featureData} and
  data values from experimental measurements stored in the \code{assayData}.}
\section{Objects from the Class}{
  Objects of the class are generally created in two different ways,
  depending on how the experimental data are provided. The first way is to
  pass an instance of \code{\link[Biobase:AssayData-class]{AssayData}} as the
  \code{dataVals} parameter in
  \code{new("trackSet", assayData,
    phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
    featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
    experimentData = new("MIAME"),
    annotation = character(), genome = "hg18")}.
  The \code{AssayData} object must contain a matrix or data.frame named
  \code{dataVals}, which should contain the experimental
  measurements. Alternatively, the experimental measurements may be
  passed as the \code{dataVals} parameter in
  \code{new("trackSet",
    phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
    featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
    experimentData = new("MIAME"),
    annotation = character(), dataVals = matrix(), genome = "hg18")}.
}
\section{Slots}{
  \describe{
    \item{\code{genome}:}{Object of class \code{"character"} identifying
      the genome to which this track pertains. Should be specified
      according to the conventions of the UCSC browser (e.g. "hg18"),
      if possible.}
    \item{\code{assayData}:}{Object of class \code{"AssayData"} holding
      the experimental measurements. }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}
      holding the experimental design matrix. }
    \item{\code{featureData}:}{Object of class
      \code{"AnnotatedDataFrame"} holding the feature information,
      generally including columns \code{chrom} (chromosome name),
      \code{start} (numeric start position),
      \code{end} (numeric end position)
      and \code{strand} (DNA strand: "+", "-" or \code{NA}). Note that
      the intervals are closed, i.e. the selected region is [start,stop].}
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} with
      experimental metadata. }
    \item{\code{annotation}:}{Object of class \code{"character"}
      identifying the annotation resource. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}
      specifying the version of the class. }
  }
}
\section{Extends}{
Class \code{"\link[Biobase:eSet-class]{eSet}"}, directly.
Class \code{"\link[Biobase:VersionedBiobase-class]{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\link[Biobase:Versioned-class]{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{\code{object[i, j]}}{
        Subsets an instance, where \code{i} is an feature index, feature
        name or \code{\linkS4class{chrid}} for selecting features and
        \code{j} is a sample index or name for selecting samples from
        the experimental data.}
    \item{\code{\link{chrid}(object)}}{Obtains an instance of
      \code{\linkS4class{chrid}}, which is a vector with
      an element representing the chromosome for each feature in this track. }
    \item{\code{\link{dataVals}(object)}}{Retrieves the
      matrix containing the experimental data values.}
    \item{\code{chrom(object)}}{Retrieves
      the chromosome identifiers of the features. }
    \item{\code{start(x)}}{Retrieves
      the start positions of the features. }
    \item{\code{end(x)}}{Retrieves the
      end positions of the features. Note that intervals are closed, so
      this should be the index of the last base to include in the feature.}
    \item{\code{strand(object)}}{Retrieves
      the strand ("+", "-" or \code{NA}) of the features. }
    \item{\code{\link{genomeSegment}(object)}}{
      obtain the genome segment spanned by this track (including gaps). }
    \item{\code{\link{trackData}(object)}}{Obtain a
      representation of this track as a \code{data.frame}. The result is
      the column-wise combination of the \code{featureData} slot with the
      \code{dataVals} matrix. The function also adds a \code{seqMid}
      column, which holds the mid-point of the segment of each
      feature. This is useful for e.g. plotting.
    }
    \item{as(trackSet, "data.frame")}{Equivalent to \code{trackData} above.}
    \item{\code{\link{trackSets}(object)}}{
      Create a \code{\linkS4class{trackSets}} instance initially containing only
      this track. }
    \item{\code{\link{export.bed}(object, con, wig = FALSE)}}{
      Export track in Browser Extended Display (\acronym{BED}) format. }
    \item{\code{\link{export.gff}(object, con,
        version = c("1", "2", "3"), source = "rtracklayer")}}{
      Export track in General Feature Format (\acronym{GFF}). }
    \item{\code{\link{export.ucsc}(object, con, subformat = c("gff", "wig"),
        name = deparse(substitute(object)), ...)}}{
      Export track in the meta format of the (\acronym{UCSC}) browser. }
    \item{\code{\link{export.wig}(object, con,
        dataFormat = c("bed", "variableStep", "fixedStep"))}}{
      export track in Wiggle (\acronym{WIG}) format. }
  }
}
\author{ Michael Lawrence }
\seealso{
  \code{\link{layTrack}} for loading a track into a genome browser.
  \code{\link{import}} to create a \code{trackSet} from a file,
  \code{\link{export}} to write a \code{trackSet} to a file.
}
\keyword{classes}
