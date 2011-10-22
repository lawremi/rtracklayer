\name{BigWigFile-class}
\docType{class}

%% Classes:
\alias{class:BigWigFile}
\alias{BigWigFile-class}

%% Constructor:
\alias{BigWigFile}

%% Accessors:
\alias{seqinfo,BigWigFile-method}

%% Import:
\alias{import.bw,BigWigFile-method}

%% Utilites:
\alias{summary,BigWigFile-method}

\title{BigWigFile objects}

\description{
  A \code{BigWigFile} object is a reference to a BigWig file. It exists
  to support methods with behavior particular to BigWig files.
}

\section{Accessor Methods}{
  In the code snippets below, \code{x} represents a \code{BigWigFile}
  object.
  
  \describe{
    \item{}{
      \code{seqinfo(x)}:
      Gets the \code{\link[GenomicRanges:Seqinfo-class]{Seqinfo}} object
      indicating the lengths of the sequences for the intervals in the
      file. No circularity or genome information is available.
    }
  }
}

\section{Import}{
  \describe{
    \item{}{
      \code{import.bw(con, selection = BigWigSelection(ranges, ...),
        ranges = con, ...)}: Imports the intervals from a big wig file
        \code{con}, according to \code{selection}, a
        \linkS4class{RangedSelection} object indicating the intervals to
        retrieve from a bigWig file. Supported types of \code{con}
        include a \code{BigWigFile} and a file name. Note that this
        retrieval is very efficient, due to the indexing of the bigWig
        format.
    }
  }
}

\section{Utilities}{
  \describe{
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
  \code{\link{import.bw}} and \code{\link{export.bw}} for reading and
  writing BigWig files, respectively.
}

\examples{
bwf <- BigWigFile(system.file("tests", "test.bw", package = "rtracklayer"))
seqinfo(bwf)
track <- import.bw(bwf, asRangedData = FALSE)
summary(bwf) # for each sequence, average all values into one
summary(bwf, range(head(track))) # just average the first few features
summary(bwf, size = GenomicRanges::seqlengths(bwf) / 10) # 10X reduction
summary(bwf, type = "min") # min instead of mean
}

\keyword{methods}
\keyword{classes}
