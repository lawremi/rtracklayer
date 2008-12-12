\name{export-tracks}
\alias{export.gff}
\alias{export.gff,RangedData-method}
\alias{export.gff1}
\alias{export.gff2}
\alias{export.gff3}
\alias{export.gff1,ANY-method}
\alias{export.gff2,ANY-method}
\alias{export.gff3,ANY-method}
\alias{export.bed}
\alias{export.bed,RangedData-method}
\alias{export.wig}
\alias{export.wig,RangedData-method}
\alias{export.ucsc}
\alias{export.ucsc,RangedData-method}
\alias{export.ucsc,RangedDataList-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Export tracks}
\description{
  These functions output \code{\link[IRanges:RangedData-class]{RangedData}} instances in
  various formats.
}
\usage{
export.gff(object, con, version = c("1", "2", "3"), source =
           "rtracklayer")
export.gff1(object, con, ...)
export.gff2(object, con, ...)
export.gff3(object, con, ...)
export.bed(object, con, wig = FALSE, color = NULL, ...)
export.wig(object, con,
           dataFormat = c("auto", "bed", "variableStep", "fixedStep"), ...)
export.ucsc(object, con, subformat = c("auto", "gff1", "wig", "bed"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ The object to export, such as a
    \code{\link[IRanges:RangedData-class]{RangedData}}. If a
    \code{\linkS4class{UCSCData}}, the track line information is
    output. In the case of \code{export.ucsc}, a
    \code{\linkS4class{RangedDataList}} 
    object with possibly multiple tracks is supported.}
  \item{con}{ The connection to which the object is exported. }
  \item{version}{ The \acronym{GFF} version, either "1", "2" or "3"
    (default is "1"). }
  \item{source}{ The source of the GFF information, for \acronym{GFF}. }
  \item{wig}{ Whether to output the \acronym{WIG} variant of
    \acronym{BED} lines, not to be used directly. }
  \item{color}{Recycled vector of colors, as interpreted by
    \code{\link{col2rgb}} for BED features. If \code{NULL}, the
    \code{color} column in the \code{featureData} is used, if any.}
  \item{dataFormat}{ The format of the data lines for \acronym{WIG}
    tracks, see references. The "auto" format uses the most efficient
    format possible.}
  \item{subformat}{ The format of the tracks within the \acronym{UCSC}
    container. If "auto", "wig" is used for numeric data, else "bed". }
  \item{\dots}{For \code{export.gff1}, \code{export.gff2} and
    \code{export.gff3}: arguments to pass to \code{export.gff}. For
    \code{export.bed} and \code{export.wig} : arguments to pass to methods. For
    \code{export.ucsc}: arguments to pass to \code{export.subformat} or
    to set on the slots of the \code{\linkS4class{TrackLine}} subclass
    corresponding to \code{subformat}.}
}
\value{
  If \code{con} is missing, a character vector containing the string
  output, otherwise nothing.
}
\references{
  \describe{
    \item{GFF1 and GFF2}{
      \url{http://www.sanger.ac.uk/Software/formats/GFF}
    }
    \item{GFF3}{\url{http://www.sequenceontology.org/gff3.shtml}}
    \item{BED}{\url{http://genome.ucsc.edu/goldenPath/help/customTrack.html\#BED}}
    \item{WIG}{\url{http://genome.ucsc.edu/goldenPath/help/wiggle.html}}
    \item{UCSC}{\url{http://genome.ucsc.edu/goldenPath/help/customTrack.html}}
  }
}
\author{ Michael Lawrence }
\seealso{ See \code{\link{export}} for the high-level interface to these
  functions. }
\examples{
  dummy <- file() # dummy file connection for demo
  track <- import(system.file("tests", "bed.wig", package = "rtracklayer"))
  ## output a track as GFF2
  export.gff(track, dummy, version = "2")
  ## equivalently
  export.gff2(track, dummy)
  ## output as WIG string in variableStep format
  wig <- export.wig(track, dummy, dataFormat = "variableStep")
  ## output multiple tracks in UCSC meta-format
  track2 <- import(system.file("tests", "v1.gff", package = "rtracklayer"))
  ## output to WIG with BED line format
  export.ucsc(RangedDataList(track, track2), dummy, subformat = "wig", dataFormat = "bed") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
