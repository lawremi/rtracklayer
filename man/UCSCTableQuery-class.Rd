\name{UCSCTableQuery-class}
\docType{class}

\alias{UCSCTableQuery-class}

% Accessors:
\alias{browserSession,UCSCTableQuery-method}
\alias{browserSession<-}
\alias{browserSession<-,UCSCTableQuery,UCSCSession-method}
\alias{trackName}
\alias{trackName,UCSCTableQuery-method}
\alias{trackName<-}
\alias{trackName<-,UCSCTableQuery-method}
\alias{tableName}
\alias{tableName,UCSCTableQuery-method}
\alias{tableName<-}
\alias{tableName<-,UCSCTableQuery-method}
\alias{range,UCSCTableQuery-method}
\alias{range<-,UCSCTableQuery-method}
\alias{names,UCSCTableQuery-method}
\alias{names<-,UCSCTableQuery-method}
\alias{trackNames,UCSCTableQuery-method}
\alias{ucscSchema}
\alias{ucscSchema,UCSCTableQuery-method}

% Query execution
\alias{tableNames}
\alias{tableNames,UCSCTableQuery-method}
\alias{getTable}
\alias{getTable,UCSCTableQuery-method}
\alias{track}
\alias{track,UCSCTableQuery-method}

% Constructor:
\alias{ucscTableQuery}
\alias{ucscTableQuery,UCSCSession-method}

% Show:
\alias{show,UCSCTableQuery-method}

\title{Querying UCSC Tables}
\description{The UCSC genome browser is backed by a large database,
  which is exposed by the Table Browser web interface. Tracks are
  stored as tables, so this is also the mechanism for retrieving tracks. The
  \code{UCSCTableQuery} class represents a query against the Table
  Browser. Storing the query fields in a formal class facilitates
  incremental construction and adjustment of a query.} 

\details{
  There are five supported fields for a table query:
  \describe{
    \item{session}{The \code{\linkS4class{UCSCSession}} instance from
      the tables are retrieved. Although all sessions are based on the
      same database, the set of user-uploaded tracks, which are represented
      as tables, is not the same, in general.
    }
    \item{trackName}{The name of a track from which to retrieve a
      table. Each track can have multiple tables. Many times there is a
      primary table that is used to display the track, while the other
      tables are supplemental. Sometimes, tracks are displayed by aggregating
      multiple tables.
    }
    \item{tableName}{The name of the specific table to retrieve. May be
      \code{NULL}, in which case the behavior depends on how the query
      is executed, see below.
    }
    \item{range}{A genome identifier, a
      \code{\link[GenomicRanges:GRanges-class]{GRanges}} or 
      a \code{\link[IRanges:RangesList-class]{RangesList}} indicating
      the portion of the table to retrieve, in genome coordinates.
      Simply specifying the genome string is the easiest way to download
      data for the entire genome, and \code{\link{GRangesForUCSCGenome}}
      facilitates downloading data for e.g. an entire chromosome.
    }
    \item{names}{Names/accessions of the desired features}
  } 

  A common workflow for querying the UCSC database is to create an
  instance of \code{UCSCTableQuery} using the \code{ucscTableQuery}
  constructor, invoke \code{tableNames} to list the available tables for
  a track, and finally to retrieve the desired table either as a
  \code{data.frame} via \code{getTable} or as a \code{RangedData} track
  via \code{track}. See the examples.

  The reason for a formal query class is to facilitate multiple queries
  when the differences between the queries are small. For example, one
  might want to query multiple tables within the track and/or same
  genomic region, or query the same table for multiple regions. The
  \code{UCSCTableQuery} instance can be incrementally adjusted for each
  new query. Some caching is also performed, which enhances performance. 
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{ucscTableQuery(x, track, range = genome(x), table = NULL,
        names = NULL)}: Creates a \code{UCSCTableQuery} with the
        \code{UCSCSession} given as \code{x} and the track name given by
        the single string \code{track}. \code{range} should be a genome
        string identifier, a \code{GRanges} instance or
        \code{RangesList} instance, and it effectively defaults to
        \code{genome(x)}. If the genome is missing, it is taken from the
        session. The table name is given by
        \code{table}, which may be a single string or
        \code{NULL}. Feature names, such as gene identifiers, may be
        passed via \code{names} as a character vector.
    }
  }
}

\section{Executing Queries}{
  Below, \code{object} is a \code{UCSCTableQuery} instance.
  
  \describe{
    \item{}{
      \code{track(object, asRangedData = TRUE)}:
      Retrieves the indicated table as a track, i.e. a \code{RangedData}
      instance. Note that not all tables are available as tracks. Pass
      \code{asRangedData = FALSE} to obtain a \code{GRanges} object.
    }
    \item{}{
      \code{getTable(object)}: Retrieves the indicated table as a
      \code{data.frame}. Note that not all tables are output in
      parseable form.
    }
    \item{}{
      \code{tableNames(object)}: Gets the names of the tables available
      for the session, track and range specified by the query.
    }
  }
}

\section{Accessor methods}{
  In the code snippets below, \code{x}/\code{object} is a
  \code{UCSCTableQuery} object.

  \describe{
    \item{}{\code{browserSession(object)},
      \code{browserSession(object) <- value}:
      Get or set the \code{UCSCSession} to query.
    }
    \item{}{\code{trackName(x)}, \code{trackName(x) <- value}: Get or
      set the single string indicating the track containing the table of
      interest.
    }
    \item{}{\code{trackNames(x)}}{List the names of the tracks available
        for retrieval for the assigned genome.}
    \item{}{\code{tableName(x)}, \code{tableName(x) <- value}: Get or
      set the single string indicating the name of the table to
      retrieve. May be \code{NULL}, in which case the table is
      automatically determined.
    }
    \item{}{\code{range(x)}, \code{range(x) <- value}: Get or set the
      \code{GRanges} indicating the portion of the table to retrieve in
      genomic coordinates. Any missing information, such as the genome
      identifier, is filled in using \code{range(browserSession(x))}. It
      is also possible to set the genome identifier string or
      a \code{RangesList}.
    }
    \item{}{\code{names(x)}, \code{names(x) <- value}: Get or set the
      names of the features to retrieve. If \code{NULL}, this filter is
      disabled.
    }
    \item{}{\code{ucscSchema(x)}: Get
      the \code{\linkS4class{UCSCSchema}} object describing the selected table.}
  }
}

\author{ Michael Lawrence }

\examples{
\dontrun{
session <- browserSession()
genome(session) <- "mm9"
trackNames(session) ## list the track names
## choose the Conservation track for a portion of mm9 chr1
query <- ucscTableQuery(session, "Conservation",
                        GRangesForUCSCGenome("mm9", "chr12",
                                             IRanges(57795963, 57815592))))
## list the table names
tableNames(query)
## get the phastCons30way track
tableName(query) <- "phastCons30way"
## retrieve the track data
track(query)
## get a data.frame summarizing the multiple alignment
tableName(query) <- "multiz30waySummary"
getTable(query)

genome(session) <- "hg18"
query <- ucscTableQuery(session, "snp129",
                        names = c("rs10003974", "rs10087355", "rs10075230"))
ucscSchema(query)
getTable(query)
}
}

\keyword{methods}
\keyword{classes}
