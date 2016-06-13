\name{BamFile-methods}

\alias{export,ANY,BamFile,ANY-method}
\alias{export,GAlignments,BamFile,ANY-method}
\alias{export,GAlignmentPairs,BamFile,ANY-method}
\alias{import,BamFile,ANY,ANY-method}

\title{Export to BAM Files}
\description{
  Methods for import and export of
  \code{\link[GenomicAlignments]{GAlignments}} or
  \code{\link[GenomicAlignments]{GAlignmentPairs}} objects from and to BAM
  files, represented as \code{\link[Rsamtools]{BamFile}} objects.
}
\usage{
\S4method{import}{BamFile,ANY,ANY}(con, format, text, paired = FALSE,
                                   use.names = FALSE,
                                   param = ScanBamParam(...),
                                   genome = NA_character_, ...)
\S4method{export}{ANY,BamFile,ANY}(object, con, format, ...)
}
\arguments{
  \item{object}{The object to export, such as
    a \code{GAlignments} or \code{GAlignmentPairs}.
  }
  \item{con}{A path, URL, connection or \code{BamFile} object.}
  \item{format}{If not missing, should be \dQuote{bam}.}
  \item{text}{Not supported.}
  \item{paired}{If \code{TRUE}, return a GAlignmentPairs object,
    otherwise a GAlignments.
  }
  \item{use.names}{Whether to parse QNAME as the names on the result.}
  \item{param}{The \code{\link[Rsamtools]{ScanBamParam}} object
    governing the import.
  }
  \item{genome}{Single string or \code{Seqinfo} object identifying the genome}
  \item{\dots}{Arguments that are passed to \code{ScanBamParam}
    if \code{param} is missing.
  }
}
\details{
  BAM fields not formally present in the \code{GAlignments[Pairs]} object
  are extracted from the metadata columns, if present; otherwise, the
  missing value, \dQuote{"."}, is output. The file is sorted and
  indexed. This can be useful for subsetting BAM files, although
  \code{\link[Rsamtools]{filterBam}} may eventually become flexible
  enough to be the favored alternative.
}
\seealso{
  The \code{\link[GenomicAlignments]{readGAlignments}} and
  \code{\link[GenomicAlignments]{readGAlignmentPairs}} functions
  for reading BAM files.
}
\examples{
     library(Rsamtools)
     ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
     gal <- import(ex1_file, param=ScanBamParam(what="flag"))
     gal.minus <- gal[strand(gal) == "-"]
\dontrun{
     export(gal, BamFile("ex1-minus.bam"))
}
}
\author{Michael Lawrence}
