\name{wigToBigWig}
\alias{wigToBigWig}
\title{
  Convert WIG to BigWig
}
\description{
  This function calls the Kent C library to efficiently convert a WIG
  file to a BigWig file, without loading the entire file into
  memory. This solves the problem where simple tools write out text WIG
  files, instead of more efficiently accessed binary BigWig files.
}
\usage{
wigToBigWig(x, seqinfo,
            dest = paste(file_path_sans_ext(x, TRUE), "bw", sep = "."))
}
\arguments{
  \item{x}{
    Path or URL to the WIG file. Connections are not supported.
  }
  \item{seqinfo}{
    \code{\linkS4class{Seqinfo}} object, describing the genome of the
    data. All BigWig files must have this defined.
  }
  \item{dest}{
    The path to which to write the BigWig file. Defaults to \code{x}
    with the extension changed to \dQuote{bw}.
  }
}

\author{
  Michael Lawrence
}

\seealso{
  \code{\link[=BigWigFile]{BigWig}} import and export support
}
