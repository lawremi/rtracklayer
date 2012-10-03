\name{BigWigFile-class}
\docType{class}

%% Classes:
\alias{class:BigWigFile}
\alias{BigWigFile-class}
\alias{class:BWFile}
\alias{BWFile-class}

%% Constructor:
\alias{BigWigFile}
\alias{BWFile}

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
\alias{export,RangedData,BigWigFile,ANY-method}
\alias{export,RleList,BigWigFile,ANY-method}

%% Utilites:
\alias{summary,BigWigFile-method}

\title{BigWig Import and Export}

\description{
  These functions support the import and export of the UCSC BigWig
  format, a compressed, binary form of WIG/BEDGraph with a spatial index
  and precomputed summaries. These functions do not work on Windows.
}

\usage{
\S4method{import}{BigWigFile,ANY,ANY}(con, format, text,
                   selection = BigWigSelection(which, ...),
                   which = con, asRangedData = TRUE, ...)
import.bw(con, ...)

\S4method{export}{ANY,BigWigFile,ANY}(object, con, format, ...)
\S4method{export}{RangedData,BigWigFile,ANY}(object, con, format,
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
    \code{RangedData} or something coercible to a \code{RangedData},
    like a \code{GRanges}.
  }
  \item{format}{If not missing, should be \dQuote{bigWig} or \dQuote{bw}
    (case insensitive).
  }
  \item{text}{Not supported.}
  \item{asRangedData}{If \code{FALSE}, a \code{GRanges} is returned,
    instead of a \code{RangedData}.
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
  A \code{RangedData} (or \code{GRanges} if the \code{asRangedData} is
  \code{FALSE}), with the score values in the \code{score} column, which
  is accessible via the \code{score} function.
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
      Gets the \code{\link[GenomicRanges:Seqinfo-class]{Seqinfo}} object
      indicating the lengths of the sequences for the intervals in the
      file. No circularity or genome information is available.
    }
    \item{}{
      \code{summary(ranges = as(seqinfo(object), "GenomicRanges"), size
        = 1L, type = c("mean", "min", "max", "coverage", "sd"),
        defaultValue = NA_real_)}: Aggregates the intervals in the file
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
    }
  }

}

\author{Michael Lawrence}

\seealso{
  \code{\link{wigToBigWig}} for converting a WIG file to BigWig.
}

\examples{
if (.Platform$OS.type != "windows") {
  test_path <- system.file("tests", package = "rtracklayer")
  test_bw <- file.path(test_path, "test.bw")

  test <- import(test_bw)
  test

  which <- GRanges(c("chr2", "chr2"), IRanges(c(1, 300), c(400, 1000)))
  import(test_bw, which = which)

\dontrun{
  test_bw_out <- file.path(tempdir(), "test_out.bw")
  export(test, test_bw_out)
}

  bwf <- BigWigFile(test_bw)
  track <- import(bwf, asRangedData = FALSE)

  seqinfo(bwf)  

  summary(bwf) # for each sequence, average all values into one
  summary(bwf, range(head(track))) # just average the first few features
  summary(bwf, size = GenomicRanges::seqlengths(bwf) / 10) # 10X reduction
  summary(bwf, type = "min") # min instead of mean
}
}

\keyword{methods}
\keyword{classes}
