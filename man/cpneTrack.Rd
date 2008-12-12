\name{cpneTrack}
\alias{cpneTrack}
\docType{data}
\title{CPNE1 SNP track}
\description{
  A \code{RangedData} object (created by the \code{GGtools} package)
  with features from a subset of the SNPs on chromosome 20 from 60
  HapMap founders in the CEU cohort. Each SNP has an associated data
  value indicating its association with the expression of the CPNE1 gene
  according to a Cochran-Armitage 1df test. The top 5000 scoring SNPs
  were selected for the track.
}
\usage{data(cpneTrack)}
\format{
  Each feature (row) is a SNP. The association test scores are
  accessible via \code{\link[IRanges]{score}}.
}
\source{
  Vince Carey and the \code{GGtools} package.
}
\examples{
data(cpneTrack)
plot(start(cpneTrack), score(cpneTrack))
}
\keyword{datasets}
