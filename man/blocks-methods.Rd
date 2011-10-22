\name{blocks-methods}
\docType{methods}
\alias{blocks}
\alias{blocks-methods}
\alias{blocks,RangedData-method}
\alias{blocks,GenomicRanges-method}
\title{Get blocks/exons}
\description{
 Obtains the block ranges (subranges, usually exons) from an object, such as a
 \code{\link[IRanges:RangedData-class]{RangedData}} imported from a BED file.
}
\usage{
blocks(x, ...)
}
\arguments{
  \item{x}{The instance from which to obtain the block/exon
    information. Currently must be a \code{RangedData} or
    \code{GenomicRanges}, presumably
    imported with \code{\link{import.bed}} or formatted with
    \code{\link{asBED}}.
  }
  \item{...}{Additional arguments for methods}
}
\details{
  For the \code{RangedData} method, there must be two columns in
  \code{x}: \code{blockStarts} and \code{blockSizes}, each field of
  which should be a comma-separated list of block starts and widths,
  respectively. This comes from the BED specification.
}
\author{Michael Lawrence}
\seealso{\code{\link{import.bed}} for importing a track from BED, which
  can store block information; \code{\link{asBED}} for coercing an
  interval dataset into a BED-like structure that can be passed to this
  function.}
\keyword{methods}
\keyword{manip}
