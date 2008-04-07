\name{export}
\alias{export}
\alias{export,ANY,ANY,character-method}
\alias{export,ANY,missing,character-method}
\alias{export,ANY,character,missing-method}
\alias{export,ANY,character,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Export objects to connections }
\description{
  Exports (serializes) an object in a given format to a given connection.
}
\usage{
export(object, con, format, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ The object to export. }
  \item{con}{ The connection to which the object is exported. If this is
  a character vector, it is assumed to be a filename and a corresponding
  file connection is created and then closed after exporting the
  object. If missing, the function will return the output as a character
  vector, rather than writing to a connection.}
  \item{format}{ The format of the output. If missing and \code{con} is
    a filename, the format is derived from the file extension. }
  \item{\dots}{ Parameters to pass to the format-specific
    export routine. }
}
\details{
  This function delegates to another function, depending on the
  specified format. The name of the delegate is of the form
  \code{export.format} where \code{format} is specified by the
  \code{format} argument.
}
\value{
  If \code{con} is missing, a character vector containing the string
  output. Otherwise, nothing is returned.
}
\author{ Michael Lawrence }
\seealso{ \code{\link{import}} for the reverse }
\examples{
  track <- import(system.file("tests", "v1.gff", package = "rtracklayer"))
  \dontrun{export(track, "my.gff", version = "3")}
  ## equivalently,
  \dontrun{export(track, "my.gff3")}
  ## or
  \dontrun{
  con <- file("my.gff3")
  export(track, con, "gff3")
  close(con)
  }
  ## or as a string
  export(track, format = "gff3")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
