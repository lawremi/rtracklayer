\name{import.gff}
\alias{import.gff}
\alias{import.gff,characterORconnection-method}
\alias{import.gff1,ANY-method}
\alias{import.gff2,ANY-method}
\alias{import.gff3,ANY-method}
\alias{import.wig,ANY-method}
\alias{import.gff1}
\alias{import.gff2}
\alias{import.gff3}
\alias{import.wig}
\alias{import.bed}
\alias{import.bed,connection-method}
\alias{import.bed,character-method}
\alias{import.bed15}
\alias{import.bed15,ANY-method}
\alias{import.bedGraph}
\alias{import.bedGraph,ANY-method}
\alias{import.ucsc}
\alias{import.ucsc,characterORconnection-method}
\alias{import.bw}
\alias{import.bw,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Importing tracks }
\description{
  These are the functions for importing
  \code{\link[IRanges:RangedData-class]{RangedData}} 
  instances from connections or text.
}
\usage{
import.gff(con, version = c("1", "2", "3"), genome = "hg18")
import.gff1(con, ...)
import.gff2(con, ...)
import.gff3(con, ...)
import.bed(con, variant = c("base", "bedGraph", "bed15"),
           trackLine = TRUE, genome = "hg18", ...)
import.bed15(con, genome = "hg18", ...)
import.bedGraph(con, genome = "hg18", ...)
import.wig(con, genome = "hg18", ...)
import.ucsc(con,
            subformat = c("auto", "gff1", "wig", "bed", "bed15", "bedGraph"),
            drop = FALSE, ...)
## not yet supported on Windows
import.bw(con, selection = BigWigSelection(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{ The connection, filename or URL from which to receive the input. }
  \item{version}{ The version of \acronym{GFF} ("1", "2" or "3"). }
  \item{genome}{ The genome to set on the imported track. }
  \item{variant}{ Variant of BED lines, not for the user.}
  \item{trackLine}{ Whether the \acronym{BED} data has a track line (it
    normally does though track lines are not mandatory). }
  \item{subformat}{ The expected subformat of the \acronym{UCSC}
    data. If "auto", automatic detection of the subformat is attempted.}
  \item{drop}{If \code{TRUE} and there is only one track in the
    \acronym{UCSC} data, return the track instead of a list.}
  \item{selection}{A \linkS4class{RangedSelection} object indicating the
    intervals to retrieve from a bigWig file. Note that this retrieval
    is very efficient, due to the indexing of the bigWig format.}
  \item{\dots}{For \code{import.gff1}, \code{import.gff2} and
    \code{import.gff3}: arguments to pass to \code{import.gff}.
    For \code{import.ucsc}: arguments to pass on to
    \code{import.subformat}. For the others, arguments to pass to methods.}
}
\value{
  An instance of \code{\link[IRanges:RangedData-class]{RangedData}} or
  one of its subclasses, except for \code{import.ucsc}, which returns a
  \code{\link[IRanges:RangedDataList-class]{RangedDataList}} instance,
  unless there is one track and the \code{drop} parameter is \code{TRUE}.
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

  # bigWig
\dontrun{
  bw <- import(system.file("tests", "test.bw", package = "rtracklayer"),
               ranges = IRanges::RangesList(chr19 = IRanges(1, 6e7)))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
