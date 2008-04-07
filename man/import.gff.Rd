\name{import.gff}
\alias{import.gff}
\alias{import.gff,ANY-method}
\alias{import.gff1,ANY-method}
\alias{import.gff2,ANY-method}
\alias{import.gff3,ANY-method}
\alias{import.wig,ANY-method}
\alias{import.gff1}
\alias{import.gff2}
\alias{import.gff3}
\alias{import.wig}
\alias{import.bed}
\alias{import.bed,ANY-method}
\alias{import.ucsc}
\alias{import.ucsc,ANY-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Importing tracks }
\description{
  These are the functions for importing \code{\linkS4class{trackSet}}
  instances from connections or text.
}
\usage{
import.gff(con, version = c("1", "2", "3"))
import.gff1(con, ...)
import.gff2(con, ...)
import.gff3(con, ...)
import.bed(con, wig = FALSE, trackLine = !wig, ...)
import.wig(con, ...)
import.ucsc(con, subformat = c("auto", "gff", "wig", "bed"),
            drop = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{ The connection from which to receive the input. }
  \item{version}{ The version of \acronym{GFF} ("1", "2" or "3"). }
  \item{wig}{ Whether the \acronym{BED} lines are expected to be
    \acronym{WIG} formatted (not for public use). }
  \item{trackLine}{ Whether the \acronym{BED} data has a track line (it
    normally does though track lines are not mandatory). }
  \item{subformat}{ The expected subformat of the \acronym{UCSC}
    data. If "auto", automatic detection of the subformat is attempted.}
  \item{drop}{If \code{TRUE} and there is only one track in the
    \acronym{UCSC} data, return the track instead of a list.}
  \item{\dots}{For \code{import.gff1}, \code{import.gff2} and
    \code{import.gff3}: arguments to pass to \code{import.gff}. For
    \code{import.bed} and \code{import.wig}: arguments to pass to methods.
    For \code{import.ucsc}: arguments to pass on to \code{import.subformat}.}
}
\value{
  An instance of \code{\linkS4class{trackSet}} or one of its subclasses,
  except for \code{import.ucsc}, which returns a
  \code{\linkS4class{trackSets}} instance, unless there is one track and
  the \code{drop} parameter is \code{TRUE}.
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
\seealso{ \code{\link{import}} for the high-level interface to these routines. }
\examples{
  # import a GFF V2 file
  gff <- import.gff(system.file("tests", "v2.gff", package = "rtracklayer"), version = "2")
  # or
  gff <- import.gff2(system.file("tests", "v2.gff", package = "rtracklayer"))

  # import a WIG file
  wig <- import.wig(system.file("tests", "bed.wig", package = "rtracklayer"))
  # or
  wig <- import.ucsc(system.file("tests", "bed.wig", package = "rtracklayer"),
                     subformat = "wig", drop = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
