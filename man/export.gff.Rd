\name{export-tracks}
\alias{export.gff}
\alias{export.gff,ANY,ANY-method}
\alias{export.gff,RangedData,characterORconnection-method}
\alias{export.gff1}
\alias{export.gff2}
\alias{export.gff3}
\alias{export.gff1,ANY-method}
\alias{export.gff2,ANY-method}
\alias{export.gff3,ANY-method}
\alias{export.bed}
\alias{export.bed,ANY,ANY-method}
\alias{export.bed,RangedData,characterORconnection-method}
\alias{export.bed,RangedDataList,ANY-method}
\alias{export.bed15}
\alias{export.bed15,ANY-method}
\alias{export.bedGraph}
\alias{export.bedGraph,ANY-method}
\alias{export.wig}
\alias{export.wig,ANY-method}
\alias{export.ucsc}
\alias{export.ucsc,ANY,ANY-method}
\alias{export.ucsc,RangedData,ANY-method}
\alias{export.ucsc,RangedDataList,ANY-method}
\alias{export.bw}
\alias{export.bw,ANY,ANY-method}
\alias{export.bw,RangedData,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Export tracks}
\description{
  These functions output
  \code{\link[IRanges:RangedData-class]{RangedData}} instances in
  various formats.
}
\usage{
export.gff(object, con, version = c("1", "2", "3"), source =
           "rtracklayer", append = FALSE, ...)
export.gff1(object, con, ...)
export.gff2(object, con, ...)
export.gff3(object, con, ...)
export.bed(object, con, variant = c("base", "bedGraph", "bed15"),
           color = NULL, append = FALSE, ...)
export.bed15(object, con, expNames = NULL, ...)
export.bedGraph(object, con, ...)
export.wig(object, con,
           dataFormat = c("auto", "variableStep", "fixedStep"), ...)
export.ucsc(object, con, subformat = c("auto", "gff1", "wig", "bed",
           "bed15", "bedGraph"), append = FALSE, ...)
## not yet supported on Windows
export.bw(object, con,
          dataFormat = c("auto", "variableStep", "fixedStep", "bedGraph"),
          seqlengths = NULL, compress = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ The object to export, such as a
    \code{\link[IRanges:RangedData-class]{RangedData}}, or anything
    coercible to a \code{RangedData}. If a
    \code{\linkS4class{UCSCData}}, the track line information is
    output. In the case of \code{export.bed15}, \code{export.bedGraph},
    \code{export.wig}, and \code{export.ucsc}, a
    \code{\link[IRanges:RangedDataList-class]{RangedDataList}} 
    object with possibly multiple tracks is supported.}
  \item{con}{ The connection to which the object is exported. }
  \item{version}{ The \acronym{GFF} version, either "1", "2" or "3"
    (default is "1"). }
  \item{source}{ The source of the GFF information, for \acronym{GFF}. }
  \item{variant}{ Which variant of BED lines to output,
    not for the user. }
  \item{color}{Recycled vector of colors, as interpreted by
    \code{\link{col2rgb}} for BED features. If \code{NULL}, the
    \code{color} column in the \code{featureData} is used, if any.}
  \item{dataFormat}{ The format of the data lines for \acronym{WIG}
    tracks, see references. The "auto" format uses the most efficient
    format possible.}
  \item{subformat}{ The format of the tracks within the \acronym{UCSC}
    container. If "auto", the type is determined from the
    trackline. If \code{object} is not a \code{UCSCData}, this
    essentially means "wig" or "bedGraph" (depending on the density) if
    there is a numeric score, else "bed".}
  \item{expNames}{ Names of the columns in \code{object} that hold the
    experimental data. Defaults to all column names,
    unless \code{object} is a \code{\linkS4class{UCSCData}}, in which
    case the \code{expNames} field is taken from the track line, if it
    exists. }
  \item{seqlengths}{The lengths of each sequence in \code{object}. If
    \code{NULL}, the chromosome lengths are retrieved for
    the \code{genome} specified on \code{object}, if possible.
  }
  \item{append}{Logical, whether to append the output to the connection}
  \item{compress}{Logical, indicating whether to compress the bigWig output}
  \item{\dots}{For \code{export.gff1}, \code{export.gff2} and
    \code{export.gff3}: arguments to pass to \code{export.gff}. For
    \code{export.bed}: arguments to pass to methods. For
    \code{export.bed15}, \code{export.bedGraph} and \code{export.wig}:
    arguments to pass to \code{export.ucsc}. For
    \code{export.ucsc}: arguments to pass to \code{export.subformat} or
    to set on the slots of the \code{\linkS4class{TrackLine}} subclass
    corresponding to \code{subformat}.}
}
\details{
  The following is some advice for choosing a file format.
  \describe{
    \item{\acronym{GFF}}{The General Feature Format is meant to
      represent any set of genomic features, with application-specific
      columns represented as \dQuote{attributes}. There are three
      principal versions (1, 2, and 3). This is a good format
      for interoperating with other genomic tools. UCSC supports GFF1,
      but it needs to be encapsulated in the UCSC metaformat,
      i.e. \code{export.ucsc(subformat = "gff1")}.}
    \item{\acronym{BED}}{The Browser Extended Display format is for
      displaying tracks in a genome browser, in particular UCSC. There
      are many options to control the appearance of the track, see
      \code{\linkS4class{GraphTrackLine}}. To output a track line
      when \code{object} is not a \code{UCSCData},
      call \code{export.ucsc(subformat = "bed")}.}
    \item{\acronym{Bed15}}{An extension of BED with 15 columns, Bed15 is
      meant to represent data from microarray experiments. Multiple
      samples/columns are supported, and the data is displayed as a
      compact heatmap. With 15 columns per feature, this format is
      probably too verbose for e.g. ChIP-seq coverage (use multiple WIG
      tracks instead).}
    \item{\acronym{bedGraph}}{A variant of BED that represents
      experimental data more compactly than \acronym{BED} and especially
      \acronym{Bed15}, although only one sample is supported. The data is
      displayed as a bar or line graph. For dense data,
      \code{WIG} is preferred.
    }
    \item{\acronym{WIG}}{The Wiggle format is meant for storing dense
      numerical data, such as the coverage from a ChIP-seq
      experiment. The data is displayed as a bar or line graph.
    }
  }
  In summary, \acronym{BED} is usually best for displaying qualitative
  features or sparse quantiative features (like ChIP-seq peaks), while
  \acronym{WIG} is usually best for displaying dense data like coverage.

  In general, columns in the \code{RangedData} are mapped to the
  column in the track format of the same name. For example, a column
  named \dQuote{itemRgb} will be mapped to the corresponding column in
  BED-formatted output, while it is ignored for other formats. Missing
  values are mapped between \code{NA} in R and the format-specific
  missing value indicator, usually \dQuote{.}. The
  following describes how the \code{RangedData} object is mapped to each
  track format. Default values for columns are given in parentheses.
  \describe{
    \item{\acronym{GFF}}{
      Maps columns named \dQuote{source} (\dQuote{rtracklayer}),
      \dQuote{feature} (\dQuote{sequence}), \dQuote{score} (\dQuote{.}),
      \dQuote{strand} (\dQuote{.}), \dQuote{frame} (\dQuote{.}), and
      (version 1 only) \dQuote{group} (\code{seqname}). In GFF versions
      2 and 3, extra columns are mapped to attributes.
    }
    \item{\acronym{BED}}{
      Maps columns named \dQuote{name} (\dQuote{.}), \dQuote{score}
      (\dQuote{.}), \dQuote{strand} (\dQuote{.}), \dQuote{thickStart}
      (\code{start}), \dQuote{thickEnd} (\code{end}), \dQuote{itemRgb}
      (\dQuote{0,0,0}), \dQuote{blockSizes}, and
      \dQuote{blockStarts}. Note that the BED field
      \dQuote{blockCounts} is derived automatically. The intervals
      specified by \dQuote{thickStart}, \dQuote{thickEnd} and
      \dQuote{blockStarts} are 0-based, half-open as in BED. Note that
      this is different from the chromosome start/end stored in the
      \code{Ranges} object (1-based, closed). The \dQuote{itemRgb}
      column should be specified in a format understood by
      \code{\link{col2rgb}}.
    }
    \item{\acronym{Bed15}}{
      In addition to the behavior for \acronym{BED} above, encodes columns
      named by the \code{expNames} parameter into the fields
      \dQuote{expCount}, \dQuote{expIds} and \dQuote{expScores}.
    }
    \item{\acronym{bedGraph}}{
      The \dQuote{score} column is used for the quantitative values.
    }
    \item{\acronym{WIG}}{
      The \dQuote{score} column is used for the quantitative values.
    }
  }
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
  ## output to WIG
  export.ucsc(RangedDataList(track, track2), dummy, subformat = "wig") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
