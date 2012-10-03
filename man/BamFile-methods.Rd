\name{BamFile-methods}

\alias{export,ANY,BamFile,ANY-method}
\alias{export,GappedAlignments,BamFile,ANY-method}
\alias{import,BamFile,ANY,ANY-method}

\title{Export to BAM Files}
\description{
  Methods for import and export of
  \code{\link[GenomicRanges]{GappedAlignments}} objects from and to BAM
  files, represented as \code{\link[Rsamtools]{BamFile}} objects.
}
\usage{
\S4method{import}{BamFile,ANY,ANY}(con, format, text, use.names = FALSE,
                                   param = ScanBamParam(...), ...)
\S4method{export}{ANY,BamFile,ANY}(object, con, format, ...)
}
\arguments{
  \item{object}{The object to export, should be
    a \code{GappedAlignments} or one of its subclasses, like
    \code{GappedReads}.
  }
  \item{con}{A path, URL, connection or \code{BamFile} object.}
  \item{format}{If not missing, should be \dQuote{bam}.}
  \item{text}{Not supported.}
  \item{use.names}{Whether to parse QNAME as the names on the result.}
  \item{param}{The \code{\link[Rsamtools]{ScanBamParam}} object
    governing the import.
  }
  \item{\dots}{Arguments that are passed to \code{ScanBamParam}
    if \code{param} is missing.
  }
}
\details{
  BAM fields not formally present in the \code{GappedAlignments} object
  are extracted from the metadata columns, if present; otherwise, the
  missing value, \dQuote{"."}, is output. The file is sorted and
  indexed. This can be useful for subsetting BAM files, although
  \code{\link[Rsamtools]{filterBam}} may eventually become flexible
  enough to be the favored alternative.
}
\seealso{
  \code{\link[GenomicRanges]{readGappedAlignments}} for reading BAM
  files into a \code{GappedAlignments}.
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
