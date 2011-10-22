\name{targets}
\Rdversion{1.1}
\alias{targets}
\docType{data}
\title{microRNA target sites}
\description{
A data frame of human microRNA target sites retrieved from MiRBase. This is a
subset of the \code{hsTargets} data frame in the \code{microRNA}
package. See the \code{rtracklayer} vignette for more details.
}
\usage{data(targets)}
\format{
  A data frame with 2981 observations on the following 6 variables.
  \describe{
    \item{\code{name}}{The miRBase ID of the microRNA.}
    \item{\code{target}}{The Ensembl ID of the targeted transcript.}
    \item{\code{chrom}}{The name of the chromosome for target site.}
    \item{\code{start}}{Target start position.}
    \item{\code{end}}{Target stop position.}
    \item{\code{strand}}{The strand of the target site, \code{"+"}, or 
         \code{"-"}.}
  }
}
\source{
The \code{microRNA} package, dataset \code{hsTargets}. Originally
MiRBase (\url{http://microrna.sanger.ac.uk/}).
}
\examples{
data(targets)
targetTrack <- with(targets, 
    GenomicData(IRanges::IRanges(start, end), 
                strand = strand, chrom = chrom))
}
\keyword{datasets}

