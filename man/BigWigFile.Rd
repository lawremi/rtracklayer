\name{BigWigFile-class}
\docType{class}

%% Classes:
\alias{class:BigWigFile}
\alias{BigWigFile-class}
\alias{class:BWFile}
\alias{BWFile-class}
\alias{class:BigWigFileList}
\alias{BigWigFileList-class}

%% Constructor:
\alias{BigWigFile}
\alias{BWFile}
\alias{BigWigFileList}

%% Accessors:
\alias{seqinfo,BigWigFile-method}

%% Import:
\alias{import.bw}
\alias{import.bw,ANY-method}
\alias{import,BigWigFile,ANY,ANY-method}

%% Export:
\alias{export.bw}
\alias{export.bw,ANY-method}
\alias{export,ANY,BigWigFile,ANY-method}
\alias{export,GenomicRanges,BigWigFile,ANY-method}
\alias{export,List,BigWigFile,ANY-method}

%% Utilites:
\alias{summary,BigWigFile-method}
\alias{path,BigWigFileList-method}
\alias{cleanupBigWigCache}

\title{BigWig Import and Export}

\description{
  These functions support the import and export of the UCSC BigWig
  format, a compressed, binary form of WIG/BEDGraph with a spatial index
  and precomputed summaries. These functions do not work on Windows.
}

\usage{
\S4method{import}{BigWigFile,ANY,ANY}(con, format, text,
                   selection = BigWigSelection(which, ...),
                   which = con,
                   as = c("GRanges", "RleList", "NumericList"), ...)
import.bw(con, ...)

\S4method{export}{ANY,BigWigFile,ANY}(object, con, format, ...)
\S4method{export}{GenomicRanges,BigWigFile,ANY}(object, con, format,
                   dataFormat = c("auto", "variableStep", "fixedStep",
                     "bedGraph"), compress = TRUE, fixedSummaries = FALSE)
export.bw(object, con, ...)
}

\arguments{
  \item{con}{A path, URL or \code{BigWigFile} object. Connections are
    not supported. For the functions ending in \code{.bw}, the file
    format is indicated by the function name. For the \code{export}
    and \code{import} methods, the format must be indicated another
    way. If \code{con} is a path, or URL, either the file
    extension or the \code{format} argument needs to be \dQuote{bigWig}
    or \dQuote{bw}.
  }
  \item{object}{The object to export, should be an \code{RleList},
    \code{IntegerList}, \code{NumericList},
    \code{GRanges} or something coercible to a \code{GRanges}.
  }
  \item{format}{If not missing, should be \dQuote{bigWig} or \dQuote{bw}
    (case insensitive).
  }
  \item{text}{Not supported.
  }
  \item{as}{Specifies the class of the return object. Default is
    \code{GRanges}, which has one range per range in the file, and a
    score column holding the value for each range. For
    \code{NumericList}, one numeric vector is returned for each range in
    the \code{selection} argument. For \code{RleList}, there is one
    \code{Rle} per sequence, and that \code{Rle} spans the entire
    sequence.
  }
  \item{selection}{A \code{\linkS4class{BigWigSelection}} object
    indicating the ranges to load.
  }
  \item{which}{A range data structure coercible to \code{IntegerRangesList},
    like a \code{GRanges}, or a \code{BigWigFile}. Only the intervals in
    the file overlapping the given ranges are returned. By default, the
    value is the \code{BigWigFile} itself. Its \code{Seqinfo} object is
    extracted and coerced to a \code{IntegerRangesList} that represents the
    entirety of the file.
  }
  \item{dataFormat}{Probably best left to \dQuote{auto}. Exists only
    for historical reasons.
  }
  \item{compress}{If \code{TRUE}, compress the data. No reason to change this.
  }
  \item{fixedSummaries}{If \code{TRUE}, compute summaries at fixed
    resolutions corresponding to the default zoom levels in the Ensembl
    genome browser (with some extrapolation): 30X, 65X, 130X, 260X,
    450X, 648X, 950X, 1296X, 4800X, 19200X. Otherwise, the resolutions
    are dynamically determined by an algorithm that computes an initial
    summary size by initializing to 10X the size of the smallest feature
    and doubling the size as needed until the size of the summary is
    less than half that of the data (or there are no further gains). It
    then computes up to 10 more levels of summary, quadrupling the size
    each time, until the summaries start to exceed the sequence size.
  }
  \item{...}{Arguments to pass down to methods to other methods. For
    import, the flow eventually reaches the \code{BigWigFile} method on
    \code{import}.
  }
}

\value{
  A \code{GRanges} (default), \code{RleList} or \code{NumericList}.
  \code{GRanges} return ranges with non-zero score values in a \code{score}
  metadata column. The length of the \code{NumericList} is the same length
  as the \code{selection} argument (one list element per range).
  The return order in the  \code{NumericList} matches the order of the
  \code{BigWigSelection} object. 
}

\section{\code{BigWigFile} objects}{
  A \code{BigWigFile} object, an extension of
  \code{\link[BiocIO:BiocFile-class]{BiocFile}} is a reference to a BigWig file.
  To cast a path, URL or connection to a \code{BigWigFile}, pass it to the
  \code{BigWigFile} constructor.

  BigWig files are more complex than most track files, and there are a
  number of methods on \code{BigWigFile} for accessing the additional
  information:

  \describe{
    \item{}{
      \code{seqinfo(x)}:
      Gets the \code{\link[GenomeInfoDb]{Seqinfo}} object
      indicating the lengths of the sequences for the intervals in the
      file. No circularity or genome information is available.
    }
    \item{}{
      \code{summary(ranges = as(seqinfo(object), "GenomicRanges"), size
        = 1L, type = c("mean", "min", "max", "coverage", "sd"),
        defaultValue = NA_real_),
	as = c("GRangesList", "RleList", "matrix"), ...}:
        Aggregates the intervals in the file
        that fall into \code{ranges}, which should be something
        coercible to \code{GRanges}.  The aggregation essentially
        compresses each sequence to a length of \code{size}. The
        algorithm is specified by \code{type}; available algorithms
        include the mean, min, max, coverage (percent sequence covered
        by at least one feature), and standard deviation. When a window
        contains no features, \code{defaultValue} is assumed. The result
        type depends on \code{as}, and can be a GRangesList, RleList or
        matrix, where the number elements (or rows) is equal to the
        length of \code{ranges}. For \code{as="matrix"}, there must be
        one unique value of \code{size}, which is equal to the number of
        columns in the result. The \code{as="matrix"} case is the only
        one that supports a \code{size} greater than the width of the
        corresponding element in \code{ranges}, where values are
	interpolated to yield the matrix result. The
        driving use case for this is visualization of coverage when the
        screen space is small compared to the viewed portion of the
        sequence. The operation is very fast, as it leverages cached
        multi-level summaries present in every BigWig file.

        If a summary statistic is not available / cannot be computed
        for a given range a warning is thrown and the defaultValue
        \code{NA_real_} is returned.
    }
  }

  When accessing remote data, the UCSC library caches data in the
  \file{/tmp/udcCache} directory. To clean the cache, call
  \code{cleanBigWigCache(maxDays)}, where any files older than
  \code{maxDays} days old will be deleted.
}

\section{\code{BigWigFileList} objects}{
  A \code{BigWigFileList()} provides a convenient way of managing a list 
  of \code{BigWigFile} instances.
}

\author{Michael Lawrence}

\seealso{
  \code{\link{wigToBigWig}} for converting a WIG file to BigWig.
}

\examples{
if (.Platform$OS.type != "windows") {
  test_path <- system.file("tests", package = "rtracklayer")
  test_bw <- file.path(test_path, "test.bw")

  ## GRanges
  ## Returns ranges with non-zero scores.
  gr <- import(test_bw)
  gr 

  which <- GRanges(c("chr2", "chr2"), IRanges(c(1, 300), c(400, 1000)))
  import(test_bw, which = which)

  ## RleList
  ## Scores returned as an RleList is equivalent to the coverage.
  ## Best option when 'which' or 'selection' contain many small ranges.
  mini <- narrow(unlist(tile(which, 50)), 2)
  rle <- import(test_bw, which = mini, as = "RleList")
  rle 

  ## NumericList
  ## The 'which' is stored as metadata:
  track <- import(test_bw, which = which, as = "NumericList")
  metadata(track)

\dontrun{
  test_bw_out <- file.path(tempdir(), "test_out.bw")
  export(test, test_bw_out)
}

  bwf <- BigWigFile(test_bw)
  track <- import(bwf)

  seqinfo(bwf)  

  summary(bwf) # for each sequence, average all values into one
  summary(bwf, range(head(track))) # just average the first few features
  summary(bwf, size = seqlengths(bwf) / 10) # 10X reduction
  summary(bwf, type = "min") # min instead of mean
  summary(bwf, track, size = 10, as = "matrix") # each feature 10 windows
}
}

\keyword{methods}
\keyword{classes}
