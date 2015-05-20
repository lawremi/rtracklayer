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

\title{BigWig Import and Export}

\description{
  These functions support the import and export of the UCSC BigWig
  format, a compressed, binary form of WIG/BEDGraph with a spatial index
  and precomputed summaries. These functions do not work on Windows.
}

\usage{
\S4method{import}{BigWigFile,ANY,ANY}(con, format, text,
                   selection = BigWigSelection(which, ...),
                   which = con, asRangedData = FALSE, asRle = FALSE, 
                   as = c("GRanges", "RleList", "NumericList"), ...)
import.bw(con, ...)

\S4method{export}{ANY,BigWigFile,ANY}(object, con, format, ...)
\S4method{export}{GenomicRanges,BigWigFile,ANY}(object, con, format,
                   dataFormat = c("auto", "variableStep", "fixedStep",
                     "bedGraph"), compress = TRUE)
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
    sequence. As of version 1.23.11 \code{import} no longer returns
    \code{RangedData}.
  }
  \item{asRangedData}{Defunct.}
  \item{asRle}{Deprecated. Use \code{as} instead.
    If \code{TRUE}, the BigWig file is assumed to contain
    contiguous ranges that define a run-length encoding of a
    vector (like coverage), and a \code{RleList} is returned.
  }
  \item{selection}{A \code{\linkS4class{BigWigSelection}} object
    indicating the ranges to load.
  }
  \item{which}{A range data structure coercible to \code{RangesList},
    like a \code{GRanges}, or a \code{BigWigFile}. Only the intervals in
    the file overlapping the given ranges are returned. By default, the
    value is the \code{BigWigFile} itself. Its \code{Seqinfo} object is
    extracted and coerced to a \code{RangesList} that represents the
    entirety of the file.
  }
  \item{dataFormat}{Probably best left to \dQuote{auto}. Exists only
    for historical reasons.
  }
  \item{compress}{If \code{TRUE}, compress the data. No reason to change this.
  } 
  \item{...}{Arguments to pass down to methods to other methods. For
    import, the flow eventually reaches the \code{BigWigFile} method on
    \code{import}. For export, the \code{RangedData}, \code{BigWigFile}
    method on \code{export} is the sink. 
  }
}

\value{
  A \code{GRanges} (default), \code{RangedData}, \code{RleList} or 
  \code{NumericList}. \code{GRanges} and \code{RangedData} return ranges
  with non-zero score values in a \code{score} metadata column. The 
  length of the \code{NumericList} is the same length as the \code{selection} 
  argument (one list element per range). The return order in the 
  \code{NumericList} matches the order of the \code{BigWigSelection} object. 
}

\section{\code{BigWigFile} objects}{
  A \code{BigWigFile} object, an extension of
  \code{\linkS4class{RTLFile}} is a reference to a BigWig file. To cast
  a path, URL or connection to a \code{BigWigFile}, pass it to the
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
        defaultValue = NA_real_), ...}: Aggregates the intervals in the file
        that fall into \code{ranges}, which should be something
        coercible to \code{GRanges}.  The aggregation essentially
        compresses each sequence to a length of \code{size}. The
        algorithm is specified by \code{type}; available algorithms
        include the mean, min, max, coverage (percent sequence covered
        by at least one feature), and standard deviation. When a window
        contains no features, \code{defaultValue} is assumed. The result
        is an \code{\link[IRanges:RleList-class]{RleList}}, with an
        element for each element in \code{ranges}. The
        driving use case for this is visualization of coverage when the
        screen space is small compared to the viewed portion of the
        sequence. The operation is very fast, as it leverages cached
        multi-level summaries present in every BigWig file.

        If a summary statistic is not available / cannot be computed
        for a given range a warning is thrown and the defaultValue
        \code{NA_real_} is returned.
    }
  }
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
}
}

\keyword{methods}
\keyword{classes}
