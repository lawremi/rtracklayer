\name{TwoBitFile-class}
\docType{class}

%% Classes:
\alias{class:TwoBitFile}
\alias{TwoBitFile-class}

%% Constructor:
\alias{TwoBitFile}

%% Accessors:
\alias{seqinfo,TwoBitFile-method}

%% Import:
\alias{import.2bit}
\alias{import.2bit,TwoBitFile-method}
\alias{import.2bit,character-method}
\alias{import.2bit,connection-method}

%% Export:
\alias{export.2bit}
\alias{export.2bit,ANY,ANY-method}
\alias{export.2bit,DNAStringSet,character-method}
\alias{export.2bit,BSgenome,character-method}

\title{2bit Files}

\description{
  The \code{export.2bit} and {import.2bit} support the export and
  import, respectively, of the UCSC 2bit compressed sequence format. The
  main advantage is speed of subsequence retrieval, as it only loads the
  sequence in the requested intervals. Compared to the FA format
  supported by Rsamtools, 2bit offers the additional feature of masking
  and also has better support in Java (and thus most genome
  browsers). The supporting \code{TwoBitFile} class is a reference to a
  TwoBit file.
}

\section{Accessor Methods}{
  In the code snippets below, \code{x} represents a \code{TwoBitFile}
  object.
  
  \describe{
    \item{}{
      \code{seqinfo(x)}:
      Gets the \code{\link[GenomicRanges:Seqinfo-class]{Seqinfo}} object
      indicating the lengths of the sequences in the
      file. No circularity or genome information is available.
    }
  }
}

\section{Import}{
  \describe{
    \item{}{
      \code{import.2bit(con, which = as(seqinfo(con), "GenomicRanges"),
      ...)}: Imports sequence from 2bit file \code{con}, which can be a
      string (path or URL) or \code{TwoBitFile} object. The sequence
      retrieval is restricted to the intervals given by \code{which},
      which should be something coercible to a \code{GRanges}. The
      returned \code{DNAStringSet} contains a \code{DNAString} for every
      interval in \code{which}.
    }
  }
}

\section{Export}{
  \describe{
    \item{}{
      \code{export.2bit(object, con, which = as(seqinfo(con),
      "GenomicRanges"), ...)}: Exports \code{object} in the two bit
      format to \code{con}, a path or URL. The \code{object} should be a
      \code{\link[Biostrings:DNAStringSet-class]{DNAStringSet}} (or
      something coercible to one) or a
      \code{\link[BSgenome:BSgenome-class]{BSgenome}} object. If a
      \code{BSgenome} object, the arguments in \code{\dots} are passed
      to \code{\link[BSgenome:BSParams-class]{bsapply}} during the export
      of each sequence.
    }
  }
}

\author{Michael Lawrence}

\examples{
tbf <- TwoBitFile(system.file("tests", "test.2bit", package = "rtracklayer"))
seqinfo(tbf)
sequence <- import.2bit(tbf) # the whole file
subrange <- IRanges::resize(as(seqinfo(tbf), "GenomicRanges"), width = 50)
subsequence <- import.2bit(tbf, which = subrange)
}

\keyword{methods}
\keyword{classes}
