\name{BigBedFile-class}
\docType{class}

%% Classes:
\alias{class:BigBedFile}
\alias{BigBedFile-class}
\alias{class:BBFile}
\alias{BBFile-class}

%% Constructor:
\alias{BigBedFile}
\alias{BBFile}

%% Accessors:
\alias{seqinfo,BigBedFile-method}

%% Import:
\alias{import.bb}
\alias{import.bb,ANY-method}
\alias{import,BigBedFile,ANY,ANY-method}

%% Export:
\alias{export.bb}
\alias{export.bb,ANY-method}
\alias{export,ANY,BigBedFile,ANY-method}
\alias{export,GenomicRanges,BigBedFile,ANY-method}

\title{BigBed Import and Export}

\description{
  These functions support the import and export of the UCSC BigBed
  format, a compressed, binary form of BED with a spatial index
  and precomputed summaries. These functions do not work on Windows.
}

\usage{
\S4method{import}{BigBedFile,ANY,ANY}(con, format, text,
                   selection = BigBedSelection(which, ...),
                   which = con, ...)
import.bb(con, ...)

\S4method{export}{ANY,BigBedFile,ANY}(object, con, format, ...)
\S4method{export}{GenomicRanges,BigBedFile,ANY}(object, con, format,
                   compress = TRUE, extraIndexes = "")
export.bb(object, con, ...)
}


\arguments{
  \item{con}{A path, URL or \code{BigBedFile} object. Connections are
    not supported. For the functions ending in \code{.bb}, the file
    format is indicated by the function name. For the \code{export}
    and \code{import} methods, the format must be indicated another
    way. If \code{con} is a path, or URL, either the file
    extension or the \code{format} argument needs to be \dQuote{bigBed}
    or \dQuote{bb}.
  }
  \item{object}{The object to export, should be \code{GRanges}.
  }
  \item{format}{If not missing, should be \dQuote{bigBed} or \dQuote{bb}
    (case insensitive).
  }
  \item{text}{Not supported.
  }
  \item{selection}{A \code{\linkS4class{BigBedSelection}} object
    indicating the ranges to load.
  }
  \item{which}{A range data structure coercible to \code{IntegerRangesList},
    like a \code{GRanges}, or a \code{BigBedFile}. Only the intervals in
    the file overlapping the given ranges are returned. By default, the
    value is the \code{BigBedFile} itself. Its \code{Seqinfo} object is
    extracted and coerced to a \code{IntegerRangesList} that represents the
    entirety of the file.
  }
  \item{compress}{If \code{TRUE}, compress the data. No reason to change this.
  }
  \item{extraIndexes}{If set, make an index on each field in a comma separated list
  }
  \item{...}{Arguments to pass down to methods to other methods. For
    import, the flow eventually reaches the \code{BigBedFile} method on
    \code{import}.
  }
}

\section{\code{BigBedFile} objects}{
  A \code{BigWigFile} object, an extension of
  \code{\linkS4class{BiocFile}} is a reference to a BigBed file. To cast
  a path, URL or connection to a \code{BigBedFile}, pass it to the
  \code{BigBedFile} constructor.

  BigBed files are more complex than most track files, and there are a
  number of methods on \code{BigBedFile} for accessing the additional
  information:

  \describe{
    \item{}{
      \code{seqinfo(x)}:
      Gets the \code{\link[GenomeInfoDb]{Seqinfo}} object
      indicating the lengths of the sequences for the intervals in the
      file. No circularity or genome information is available.
    }
  }

  When accessing remote data, the UCSC library caches data in the
  \file{/tmp/udcCache} directory. To clean the cache, call
  \code{cleanBigBedCache(maxDays)}, where any files older than
  \code{maxDays} days old will be deleted.
}

\author{Michael Lawrence}

\examples{
if (.Platform$OS.type != "windows") {
  test_path <- system.file("tests", package = "rtracklayer")
  test_bb <- file.path(test_path, "test.bb")

  ## Returns ranges with all fields
  gr <- import(test_bb)
  gr

  ## Retuns ranges only for 'chr10'
  ## between 180185-180185 with all fields
  which <- GRanges(c("chr10"), IRanges(c(180185, 180185)))
  import(test_bb, which = which)

  ## Retuns ranges only for 'chr10' 
  ## between 180185-180185 with name and peak fields
  selection <- BigBedSelection(which, colnames = c("name", "peak"))
  import(test_bb, selection = selection)

\dontrun{
  test_bb_out <- file.path(tempdir(), "test_out.bb")
  export(test, test_bb_out)

  ## make an index for 'name'
  test_bb_out <- file.path(tempdir(), "test_out.bb")
  export(test, test_bb_out, extraIndexes = "name")
}
}
}

\keyword{methods}
\keyword{classes}
