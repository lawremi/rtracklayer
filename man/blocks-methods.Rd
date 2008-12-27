\name{blocks-methods}
\docType{methods}
\alias{blocks}
\alias{blocks-methods}
\alias{blocks,RangedData-method}
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
    information. Currently must be a \code{RangedData}, presumably
    imported with \code{\link{import.bed}}.
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
  can store block information.}
\keyword{methods}
\keyword{manip}
