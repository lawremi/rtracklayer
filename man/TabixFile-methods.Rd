\name{TabixFile-methods}

\alias{import,TabixFile,character,ANY-method}
\alias{import,TabixFile,missing,ANY-method}

\title{TabixFile Import}

\description{
  These methods support the import of
  \code{Rsamtools:TabixFile}{TabixFile} objects. These are generally
  useful when importing tabix-indexed files that have a non-standard
  format (i.e., not BED nor GFF). This relies on the tabix header,
  which indicates the columns in the file that correspond to the
  chromosome, start and end. The BED and GFF parsers handle tabix
  transparently. 
}

\usage{
\S4method{import}{TabixFile,character,ANY}(con, format, text,
                   which = if (is.na(genome)) NULL
                           else as(seqinfoForGenome(genome), "GenomicRanges"),
                   genome = NA, header = TRUE, ...)
\S4method{import}{TabixFile,missing,ANY}(con, format, text, ...)
}
\arguments{
  \item{con}{A \code{TabixFile} object}
  \item{format}{If any known format, like \dQuote{bed} or \dQuote{gff}
    (or one of their variants), then the appropriate parser is
    applied. If any other value, then the tabix header is consulted for
    the format. By default, this is taken from the file extension.
  }
  \item{text}{Ignored.}
  \item{which}{A range data structure coercible to \code{RangesList},
    like a \code{GRanges}. Only the intervals in the file overlapping the
    given ranges are returned. The default is to use the range over the
    entire genome given by \code{genome}, if specified.
  }
  \item{genome}{The identifier of a genome, or \code{NA} if
    unknown. Typically, this is a UCSC identifier like \dQuote{hg19}. An
    attempt will be made to derive the \code{seqinfo} on the return
    value using either an installed BSgenome package or UCSC, if network
    access is available.
  }
  \item{header}{If \code{TRUE}, then the header in the indexed file,
    which might include a track line, is sent to the parser. Otherwise,
    the initial lines are skipped, according to the \code{skip} field in
    the tabix index header.
  }
  \item{...}{Extra arguments to pass to the underlying import routine,
    which for non-standard formats is \code{\link{read.table}}.
  }
}
\value{
  A \code{RangedData} or \code{GRanges}, depending on arguments.
}
\references{
  \url{http://samtools.sourceforge.net/tabix.shtml}
}
\seealso{
  \code{\link[Rsamtools]{scanTabix}} and friends
}
\author{Michael Lawrence}
