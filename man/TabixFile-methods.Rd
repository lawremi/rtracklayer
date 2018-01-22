\name{TabixFile-methods}

\alias{import,TabixFile,ANY,ANY-method}

\alias{exportToTabix}
\alias{exportToTabix,ANY,character-method}

% undocumented/internal for now
\alias{export,GenomicRanges,TabSeparatedFile,ANY-method}
\alias{export,ANY,TabSeparatedFile,ANY-method}

\title{TabixFile Import/Export}

\description{
  These methods support the import and export of
  \code{Rsamtools:TabixFile}{TabixFile} objects. These are generally
  useful when working with tabix-indexed files that have a non-standard
  format (i.e., not BED nor GFF), as well as exporting an object with
  arbitrary columns (like a GRanges) to an indexed, tab-separated
  file. This relies on the tabix header, which indicates the columns in
  the file that correspond to the chromosome, start and end. The BED and
  GFF parsers handle tabix transparently.
}

\usage{
\S4method{import}{TabixFile,ANY,ANY}(con, format, text,
                   which = if (is.na(genome)) NULL
                           else as(seqinfoForGenome(genome), "GenomicRanges"),
                   genome = NA, header = TRUE, ...)
exportToTabix(object, con, ...)
}
\arguments{
  \item{con}{For \code{import}, a \code{TabixFile} object;
    for \code{exportToTabix}, a string naming the destination file.
  }
  \item{object}{The object to export. It is coerced to a
    \code{data.frame}, written to a tab-separated file, and indexed with
    tabix for efficient range-based retrieval of the data using
    \code{import}.
  }
  \item{format}{If any known format, like \dQuote{bed} or \dQuote{gff}
    (or one of their variants), then the appropriate parser is
    applied. If any other value, then the tabix header is consulted for
    the format. By default, this is taken from the file extension.
  }
  \item{text}{Ignored.}
  \item{which}{A range data structure coercible to \code{IntegerRangesList},
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
    which for non-standard formats is \code{\link{read.table}} or
    \code{\link{write.table}}.
  }
}
\value{
  For \code{import}, a \code{GRanges} object.

  For \code{exportToTabix}, a \code{TabixFile} object that is directly
  passable to \code{import}.
}
\references{
  \url{http://samtools.sourceforge.net/tabix.shtml}
}
\seealso{
  \code{\link[Rsamtools]{scanTabix}} and friends
}
\author{Michael Lawrence}
