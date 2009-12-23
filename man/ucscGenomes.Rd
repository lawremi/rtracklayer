\name{ucscGenomes}
\alias{ucscGenomes}
\title{Get available genomes on UCSC}
\description{
  Get a \code{data.frame} describing the available UCSC genomes.
}
\usage{
ucscGenomes()
}
\value{
  A \code{data.frame} with the following columns:
  \description{
    \item{db}{UCSC DB identifier (e.g. "hg18")}
    \item{organism}{The name of the species (e.g. "Human")}
    \item{date}{The date the genome was built}
    \item{name}{The official name of the genome build}
  }
}
\seealso{\code{\linkS4class{UCSCSession}} for details on specifying the
  genome.}
\examples{
ucscGenome()
}
\author{ Michael Lawrence }
\keyword{interface}
