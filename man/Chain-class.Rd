\name{Chain-class}
\docType{class}

%% Classes:
\alias{class:Chain}
\alias{Chain-class}
\alias{ChainBlock-class}
\alias{class:ChainBlock}
\alias{class:ChainFile}
\alias{ChainFile-class}

%% Constructor:
\alias{import.chain}
\alias{import.chain,ANY-method}
\alias{import,ChainFile,ANY,ANY-method}
\alias{ChainFile}

%% Generics and methods:
\alias{ranges,ChainBlock-method}
\alias{offset,ChainBlock-method}
\alias{score,ChainBlock-method}
\alias{space,ChainBlock-method}
\alias{reversed}
\alias{reversed,ChainBlock-method}

\title{Chain objects}

\description{
  A \code{Chain} object represents a UCSC chain alignment, typically
  imported from a \code{chain} file, and is essentially a list of
  \code{ChainBlock} objects. Each \code{ChainBlock} has a corresponding
  chromosome (its name in the list) and is a run-length
  encoded alignment, mapping a set of intervals on that chromosome to
  intervals on the same or other chromosomes. 
}

\note{
  A chain file essentially details many local alignments, so it is
  possible for the "from" ranges to map to overlapping regions in the
  other sequence. The "from" ranges are guaranteed to be disjoint (but
  do not necessarily cover the entire "from" sequence).
}

\section{Accessor Methods}{
  In the code snippets below, \code{x} and \code{object} are
  \code{ChainBlock} objects.
  
  \describe{
    \item{}{
      \code{ranges(x)}:
      Get the \code{\link[IRanges:Ranges-class]{Ranges}} object holding
      the starts and ends of the "from" ranges. Each range is a
      contiguous block of positions aligned without gaps to the other
      sequence.
    }
    \item{}{
      \code{offset(x)}: Integer offset from the "from" start to the
      "end" start (which could be in another chromosome).
    }
    \item{}{
      \code{score(x)}: The score for each mapping.
    }
    \item{}{
      \code{space(x)}: The space (chromosome) of the "to" range.
    }
    \item{}{
      \code{reversed(x)}: Whether the mapping inverts the region, i.e.,
      the alignment is between different strands.
    }
  }
}

\section{Import}{
  A \code{Chain} object can be loaded from a UCSC chain format file
  simply by passing the path \code{import} function. If the file
  extension is not \dQuote{chain}, then either pass \dQuote{chain} to
  the \code{format} argument, or cast the path to a \code{ChainFile}
  object. The \code{import.chain} function is provided as a (slight)
  convenience. It is documented below, along with the extra
  \code{exclude} argument to the \code{import} method.
  
  \describe{
    \item{}{
      \code{import.chain(con, exclude = "_", ...)}:
      Imports a chain file named \code{con} as a \code{Chain} object, a
      list of \code{ChainBlock}s. Alignments for chromosomes matching the
      \code{exclude} pattern are not imported.
    }
  }
}

\author{Michael Lawrence}

\seealso{
  \code{\link{liftOver}} for performing lift overs using a chain alignment
}

\keyword{methods}
\keyword{classes}
