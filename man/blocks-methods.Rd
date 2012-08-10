\name{blocks-methods}
\docType{methods}
\alias{blocks}
\alias{blocks-methods}
\alias{blocks,RangedData-method}
\alias{blocks,GenomicRanges-method}
\title{Get blocks/exons}
\description{
 Obtains the block ranges (subranges, usually exons) from an object,
 such as a \code{\link[IRanges:RangedData-class]{RangedData}} imported
 from a BED file.
}
\usage{
blocks(x, ...)
}
\arguments{
  \item{x}{The instance from which to obtain the block/exon
    information. Currently must be a \code{RangedData} or
    \code{GenomicRanges}, with a value column of name \dQuote{blocks}
    and of type \code{RangesList}. Such an object is returned by
    \code{\link{import.bed}} and \code{\link{asBED}}.
  }
  \item{...}{Additional arguments for methods}
}
\value{
  A \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} with an
  element for each range in \code{x}. The original block ranges are
  relative to the start of the containing range, so the returned ranges
  are shifted to absolute coordinates. The seqname and strand are
  inherited from the containing range.
}
\author{Michael Lawrence}
\seealso{\code{\link{import.bed}} for importing a track from BED, which
  can store block information; \code{\link{asBED}} for coercing a
  \code{GenomicRanges} into a BED-like structure that can be passed to
  this function.
}
\keyword{methods}
\keyword{manip}
