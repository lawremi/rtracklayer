\name{genomeSegment-class}
\docType{class}
\alias{genomeSegment-class}
\alias{merge,genomeSegment,genomeSegment-method}

\alias{genome}
\alias{genome,genomeSegment-method}
\alias{chrom,genomeSegment-method}
\alias{start,genomeSegment-method}
\alias{end,genomeSegment-method}
\alias{genome<-}
\alias{genome<-,genomeSegment-method}
\alias{chrom<-}
\alias{chrom<-,genomeSegment-method}
\alias{start<-,genomeSegment-method}
\alias{end<-,genomeSegment-method}

\title{Class "genomeSegment"}
\description{ An object representing a segment of a genome. }
\section{Objects from the Class}{
  Objects are normally created by a call of the following form:
  \code{\link{genomeSegment}(genome, chrom, start, end, segment)}. All
  parameters are optional and correspond to slots of the new
  \code{genomeSegment} instance, except for the \code{segment} parameter,
  which serves as the base for the new instance. The slots of the new
  segment match those of \code{segment}, unless the call specifies other
  parameters, which override the corresponding slots in \code{segment}.
}
\section{Slots}{
  \describe{
    \item{\code{genome}:}{Object of class \code{"character"} giving the
      name of the genome (e.g. "hg18"). }
    \item{\code{chrom}:}{Object of class \code{"character"} giving the
      chromosome name (e.g. "chr22"). }
    \item{\code{start}:}{Object of class \code{"numeric"} indicating the
      start position of the segment on the chromosome.}
    \item{\code{end}:}{Object of class \code{"numeric"} indicating the
      stop position of the segment on the chromosome. }
  }
}
\section{Methods}{
  \describe{
    \item{\code{merge(x, y)}}{
      Merges two instances. The result has the same
      slot values as \code{x}, except the non-empty slots in
      \code{y} override those from \code{x}.
    }
    \item{\code{genome(object)}}{
      Gets the genome identifier.
    }
    \item{\code{genome(object) <- value}}{
      Sets the genome identifier.
    }
    \item{\code{chrom(object)}}{
      Gets the \code{chrom}osome identifier.
    }
    \item{\code{chrom(object) <- value}}{
      Sets the \code{chrom}osome identifier.
    }
    \item{\code{start(object)}}{
      Gets the start position.
    }
    \item{\code{start(object) <- value}}{
      Sets the start position.
    }
    \item{\code{end(object)}}{
      Gets the end position.
    }
    \item{\code{end(object) <- value}}{
      Sets the end position.
    }
  }
}
\author{ Michael Lawrence }

\seealso{
  \code{\link{genomeSegment}} for obtaining instances of this class.
}
\keyword{classes}
