\name{TwoBitFile-class}
\docType{class}

%% Classes:
\alias{class:TwoBitFile}
\alias{TwoBitFile-class}
\alias{class:2BitFile}
\alias{2BitFile-class}

%% Constructor:
\alias{TwoBitFile}
\alias{2BitFile}

%% Accessors:
\alias{seqinfo,TwoBitFile-method}

%% Import:
\alias{import.2bit}
\alias{import.2bit,ANY-method}
\alias{import,TwoBitFile,ANY,ANY-method}

%% Export:
\alias{export.2bit}
\alias{export.2bit,ANY-method}
\alias{export,DNAStringSet,TwoBitFile,ANY-method}
\alias{export,DNAStringSet,character,ANY-method}
\alias{export,ANY,TwoBitFile,ANY-method}
\alias{export,BSgenome,TwoBitFile,ANY-method}
\alias{export,BSgenome,character,ANY-method}

\title{2bit Files}

\description{
  These functions support the import and export of the UCSC 2bit
  compressed sequence format. The main advantage is speed of subsequence
  retrieval, as it only loads the sequence in the requested
  intervals. Compared to the FA format supported by Rsamtools, 2bit
  offers the additional feature of masking and also has better support
  in Java (and thus most genome browsers). The supporting
  \code{TwoBitFile} class is a reference to a TwoBit file.
}

\usage{
\S4method{import}{TwoBitFile,ANY,ANY}(con, format, text,
           which = as(seqinfo(con), "GenomicRanges"), ...)
import.2bit(con, ...)

\S4method{export}{ANY,TwoBitFile,ANY}(object, con, format, ...)
\S4method{export}{DNAStringSet,TwoBitFile,ANY}(object, con, format)
\S4method{export}{DNAStringSet,character,ANY}(object, con, format, ...)
\S4method{export}{BSgenome,TwoBitFile,ANY}(object, con, format, ...)
\S4method{export}{BSgenome,character,ANY}(object, con, format, ...)
export.2bit(object, con, ...)
}

\arguments{
  \item{con}{A path, URL or \code{TwoBitFile} object. Connections are
    not supported. For the functions ending in \code{.2bit}, the file
    format is indicated by the function name. For the \code{export}
    and \code{import} methods, the format must be indicated another
    way. If \code{con} is a path, or URL, either the file
    extension or the \code{format} argument needs to be \dQuote{twoBit}
    or \dQuote{2bit}.
  }
  \item{object}{The object to export, either a \code{DNAStringSet},
    \code{BSgenome} and something coercible to a \code{DNAStringSet},
    like a character vector.
  }
  \item{format}{If not missing, should be \dQuote{twoBit} or \dQuote{2bit}
    (case insensitive).
  }
  \item{text}{Not supported.}
  \item{which}{A range data structure coercible to \code{RangesList},
    like a \code{GRanges}, or a \code{TwoBitFile}. Only the intervals in
    the file overlapping the given ranges are returned. By default, the
    value is the \code{TwoBitFile} itself. Its \code{Seqinfo} object is
    extracted and coerced to a \code{RangesList} that represents the
    entirety of the file.
  }
  \item{...}{Arguments to pass down to methods to other methods. For
    import, the flow eventually reaches the \code{TwoBitFile} method on
    \code{import}. For export, the \code{TwoBitFile} methods on
    \code{export} are the sink. The \code{BSgenome} method forwards
    arguments to \code{\link[BSgenome]{bsapply}}.
  }
}

\value{
  For import, a \code{DNAStringSet}.
}

\section{\code{TwoBitFile} objects}{
  A \code{TwoBitFile} object, an extension of
  \code{\linkS4class{RTLFile}} is a reference to a TwoBit file. To cast
  a path, URL or connection to a \code{TwoBitFile}, pass it to the
  \code{TwoBitFile} constructor.

  A TwoBit file embeds the sequence information, which can be retrieved
  with the following:
  \describe{
    \item{}{
      \code{seqinfo(x)}:
      Gets the \code{\link[GenomicRanges:Seqinfo-class]{Seqinfo}} object
      indicating the lengths of the sequences for the intervals in the
      file. No circularity or genome information is available.
    }
  }
}

\author{Michael Lawrence}

\examples{
  test_path <- system.file("tests", package = "rtracklayer")
  test_2bit <- file.path(test_path, "test.2bit")

  test <- import(test_2bit)
  test

  test_2bit_file <- TwoBitFile(test_2bit)
  import(test_2bit_file) # the whole file
  
  which_range <- IRanges(c(10, 40), c(30, 42))
  which <- GRanges(names(test), which_range)
  import(test_2bit, which = which)

  seqinfo(test_2bit_file)

\dontrun{
  test_2bit_out <- file.path(tempdir(), "test_out.2bit")
  export(test, test_2bit_out)

  ## just a character vector
  test_char <- as.character(test)
  export(test_char, test_2bit_out)
}
}

\keyword{methods}
\keyword{classes}
