\name{UCSCSchema-class}
\docType{class}

\alias{UCSCSchema-class}
\alias{genome,UCSCSchema-method}
\alias{tableName,UCSCSchema-method}
\alias{genome,UCSCSchema-method}
\alias{formatDescription}
\alias{formatDescription,UCSCSchema-method}
\alias{nrow,UCSCSchema-method}
\alias{ucscSchema,UCSCSchemaDescription-method}

% Accessors:

\title{UCSC Schema}
\description{This is a preliminary class that describes a table in the
  UCSC database. The description includes the table name, corresponding
  genome, row count, and a textual description of the format. In the
  future, we could provide more table information, like the links and
  sample data frame. This is awaiting a use-case.} 

\section{Accessor methods}{
  In the code snippets below, \code{x}/\code{object} is a
  \code{UCSCSchema} object.

  \describe{
    \item{}{\code{genome(x)}: Get the genome for the table. }
    \item{}{\code{tableName(x)}: Get the name of the table. }
    \item{}{\code{nrow(x)}: Get the number of rows in the table. }
    \item{}{\code{formatDescription(x)}: Get a textual description of
      the table format. }
  }
}

\author{ Michael Lawrence }

\examples{
\dontrun{
session <- browserSession()
genome(session) <- "mm9"
query <- ucscTableQuery(session, "knownGene")
schema <- ucscSchema(query)
nrow(schema)
}
}

\keyword{methods}
\keyword{classes}
