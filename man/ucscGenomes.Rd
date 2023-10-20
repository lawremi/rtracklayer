\name{ucscGenomes}
\alias{ucscGenomes}
\title{Get available genomes on UCSC}
\description{
  Get a \code{data.frame} describing the available UCSC genomes.
}
\usage{
ucscGenomes(organism=FALSE)
}
\arguments{
  \item{organism}{A \code{logical(1)} indicating whether scientific name
    should be appended.}
}
\value{
  A \code{data.frame} with the following columns:
  \item{db}{UCSC DB identifier (e.g. "hg18")}
  \item{species}{The name of the species (e.g. "Human")}
  \item{date}{The date the genome was built}
  \item{name}{The official name of the genome build}
  \item{organism}{The scientific name of the species (e.g. "Homo sapiens")}
}
\seealso{\code{\linkS4class{UCSCSession}} for details on specifying the
  genome.}
\details{
  For populating the organism column, the web url
  \url{https://genome.ucsc.edu/cgi-bin} is scraped for every assembly
  version to get the scientific name.
}
\examples{
ucscGenomes()
}
\author{ Michael Lawrence }

\keyword{interface}
