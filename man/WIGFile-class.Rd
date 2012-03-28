\name{WIGFile-class}
\docType{class}

%% Classes:
\alias{class:WIGFile}
\alias{WIGFile-class}

%% Constructor:
\alias{WIGFile}

%% Import:
\alias{import,WIGFile,ANY,ANY-method}
\alias{import.wig}
\alias{import.wig,ANY-method}

%% Export:
\alias{export,ANY,WIGFile,ANY-method}
\alias{export,RangedData,WIGFile,ANY-method}
\alias{export,RangedDataList,WIGFile,ANY-method}
\alias{export,UCSCData,WIGFile,ANY-method}
\alias{export.wig}
\alias{export.wig,ANY-method}

\title{WIG Import and Export}

\description{
  These functions support the import and export of the UCSC WIG (Wiggle)
  format.
}

\usage{
\S4method{import}{WIGFile,ANY,ANY}(con, format, text, genome = NA,
                    asRangedData = TRUE, trackLine = TRUE, which = NULL, ...)
import.wig(con, ...)

\S4method{export}{ANY,WIGFile,ANY}(object, con, format, ...)
\S4method{export}{RangedData,WIGFile,ANY}(object, con, format,
                   dataFormat = c("auto", "variableStep", "fixedStep"),
                   writer = .wigWriter, append = FALSE, ...)
\S4method{export}{RangedDataList,WIGFile,ANY}(object, con, format, ...)
\S4method{export}{UCSCData,WIGFile,ANY}(object, con, format,
                   trackLine = TRUE, ...)
export.wig(object, con, ...)
}

\arguments{
  \item{con}{A path, URL, connection or \code{WIGFile} object. For the
    functions ending in \code{.wig}, the file format is indicated by the
    function name. For the base \code{export} and \code{import}
    functions, the format must be indicated another way. If \code{con}
    is a path, URL or connection, either the file extension or the
    \code{format} argument needs to be \dQuote{wig}. Compressed files
    (\dQuote{gz}, \dQuote{bz2} and \dQuote{xz}) are handled
    transparently.
  }
  \item{object}{The object to export, should be a \code{RangedData} or
    something coercible to a \code{RangedData}, like a
    \code{GRanges}. For exporting multiple tracks, in the UCSC track
    line metaformat, pass a \code{RangedDataList}, or something
    coercible to one, like a \code{GenomicRangesList}.
  }
  \item{format}{If not missing, should be \dQuote{wig}.
  }
  \item{text}{If \code{con} is missing, a character vector to use as the
    input
  }
  \item{trackLine}{Whether to parse/output a UCSC track line. An
    imported track line will be stored in a \code{\linkS4class{TrackLine}}
    object, as part of the returned \code{\linkS4class{UCSCData}}.
  }
  \item{genome}{The identifier of a genome, or \code{NA} if
    unknown. Typically, this is a UCSC identifier like \dQuote{hg19}. An
    attempt will be made to derive the \code{seqinfo} on the return
    value using either an installed BSgenome package or UCSC, if network
    access is available.
  }
  \item{asRangedData}{If \code{FALSE}, a \code{GRanges} is returned,
    instead of a \code{RangedData}.
  }
  \item{which}{A range data structure like \code{RangesList}
    or \code{GRanges}. Only the intervals in the file overlapping the
    given ranges are returned. This is inefficient; use BigWig for
    efficient spatial queries.
  }
  \item{append}{If \code{TRUE}, and \code{con} points to a file path,
    the data is appended to the file. Obviously, if \code{con} is a
    connection, the data is always appended.
  }
  \item{dataFormat}{Probably best left to \dQuote{auto}. Exists only
    for historical reasons.
  }
  \item{writer}{Function for writing out the blocks; for internal use
    only.
  } 
  \item{...}{Arguments to pass down to methods to other methods. For
    import, the flow eventually reaches the \code{WIGFile} method on
    \code{import}. For export, the \code{RangedData}, \code{WIGFile}
    method on \code{export} is the sink. When \code{trackLine} is
    \code{TRUE}, the arguments are passed through \code{export.ucsc}, so
    track line parameters are supported.
  }
}

\value{
  A \code{RangedData} (or \code{GRanges} if the \code{asRangedData} is
  \code{FALSE}), with the score values in the \code{score} column, which
  is accessible via the \code{score} function.
}

\details{
  The WIG format is a text-based format for efficiently representing a
  dense genome-scale score vector. It encodes, for each feature, a range
  and score. Features from the same sequence (chromosome) are grouped
  together into a block, with a single block header line indicating the
  chromosome. There are two block formats: fixed step and variable
  step. For fixed step, the number of positions (or step) between
  intervals is the same across an entire block. For variable step, the
  start position is specified for each feature. For both fixed and
  variable step, the span (or width) is specified in the header and thus
  must be the same across all features. This requirement of uniform
  width dramatically limits the applicability of WIG. For scored
  features of variable width, consider \link[=BEDGraphFile]{BEDGraph} or
  \link[=BigWigFile]{BigWig}, which is generally preferred over both WIG
  and BEDGraph. To efficiently convert an existing WIG or BEDGraph file
  to BigWig, call \code{\link{wigToBigWig}}. Neither WIG, BEDGraph nor
  BigWig allow overlapping features.
}

\section{WIGFile objects}{
  The \code{WIGFile} class extends \code{\linkS4class{RTLFile}} and is a
  formal represention of a resource in the WIG format.
  To cast a path, URL or connection to a \code{WIGFile}, pass it to
  the \code{WIGFile} constructor.
}

\author{Michael Lawrence}

\references{
  \url{http://genome.ucsc.edu/goldenPath/help/wiggle.html}
}

\examples{
  test_path <- system.file("tests", package = "rtracklayer")
  test_wig <- file.path(test_path, "step.wig")

  ## basic import calls
  test <- import(test_wig)
  test
  import.wig(test_wig)
  test_wig_file <- WIGFile(test_wig)
  import(test_wig_file)
  test_wig_con <- file(test_wig)
  import(test_wig_con, format = "wig")
  close(test_wig_con)
  test_wig_con <- file(test_wig)
  import(WIGFile(test_wig_con))
  close(test_wig_con)

  ## various options
  import(test_wig, genome = "hg19")
  import(test_wig, trackLine = FALSE)
  import(test_wig, asRangedData = FALSE)
  which <- ranges(test[3:4,])
  import(test_wig, which = which)

\dontrun{
  ## basic export calls
  test_wig_out <- file.path(tempdir(), "test.wig")
  export(test, test_wig_out)
  export.wig(test, test_wig_out)
  test_foo_out <- file.path(tempdir(), "test.foo")
  export(test, test_foo_out, format = "wig")
  test_wig_out_file <- WIGFile(test_wig_out)
  export(test, test_wig_out_file)

  ## appending
  test2 <- initialize(test,
                      trackLine = initialize(test@trackLine,
                        name = "test2"))
  export(test2, test_wig_out_file, append = TRUE)

  ## passing track line parameters
  export(test, test_wig_out, name = "test2")

  ## no track line
  export(test, test_wig_out, trackLine = FALSE)
  
  ## gzip
  test_wig_gz <- paste(test_wig_out, ".gz", sep = "")
  export(test, test_wig_gz)
}
}

\keyword{methods}
\keyword{classes}
