\name{UCSCTableQuery-class}
\docType{class}

\alias{UCSCTableQuery-class}

% Accessors:
\alias{genome,UCSCTableQuery-method}
\alias{genome<-,UCSCTableQuery-method}
\alias{hubUrl}
\alias{hubUrl,UCSCTableQuery-method}
\alias{hubUrl<-,UCSCTableQuery-method}
\alias{tableName}
\alias{tableName,UCSCTableQuery-method}
\alias{tableName<-}
\alias{tableName<-,UCSCTableQuery-method}
\alias{range,UCSCTableQuery-method}
\alias{range<-,UCSCTableQuery-method}
\alias{names,UCSCTableQuery-method}
\alias{names<-,UCSCTableQuery-method}
\alias{ucscSchema}
\alias{ucscSchema,UCSCTableQuery-method}
\alias{ucscTables}

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
\alias{ucscTableQuery,character-method}

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
  There are six supported fields for a table query:
  \describe{
    \item{provider}{
      The provider should be a session, a genome identifier, or a TrackHub URI.
      \code{session}: The \code{\linkS4class{UCSCSession}} instance from
      the tables are retrieved. Although all sessions are based on the
      same database, the set of user-uploaded tracks, which are represented
      as tables, is not the same, in general.
    }
    \item{tableName}{The name of the specific table to retrieve. May be
      \code{NULL}, in which case the behavior depends on how the query
      is executed, see below.
    }
    \item{range}{A genome identifier, a
      \code{\link[GenomicRanges]{GRanges}} or 
      a \code{\link[IRanges]{IntegerRangesList}} indicating
      the portion of the table to retrieve, in genome coordinates.
      Simply specifying the genome string is the easiest way to download
      data for the entire genome, and \code{\link{GRangesForUCSCGenome}}
      facilitates downloading data for e.g. an entire chromosome.
    }
    \item{hubUrl}{The URI of the specific TrackHub}
    \item{genome}{A genome identifier of the specific TrackHub, only need to provide it if the provider is up of TrackHub URI.}
    \item{names}{Names/accessions of the desired features}
  } 

  A common workflow for querying the UCSC database is to create an
  instance of \code{UCSCTableQuery} using the \code{ucscTableQuery}
  constructor, invoke \code{tableNames} to list the available tables for
  a track, and finally to retrieve the desired table either as a
  \code{data.frame} via \code{getTable} or as a track
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
      \code{ucscTableQuery(x, range = seqinfo(x), table = NULL,
        names = NULL, hubUrl = NULL, genome = NULL)}: Creates a \code{UCSCTableQuery} with the
        \code{UCSCSession}, genome identifier or TrackHub URI given as \code{x} and
        the table name given by the single string \code{table}. \code{range} should
        be a genome string identifier, a \code{GRanges} instance or
        \code{IntegerRangesList} instance, and it effectively defaults to
        \code{genome(x)}. If the genome is missing, it is taken from the
        provider. Feature names, such as gene identifiers, may be
        passed via \code{names} as a character vector.
    }
  }
}

\section{Executing Queries}{
  Below, \code{object} is a \code{UCSCTableQuery} instance.
  
  \describe{
    \item{}{
      \code{track(object)}:
      Retrieves the indicated table as a track, i.e. a \code{GRanges}
      object. Note that not all tables are available as tracks.
    }
    \item{}{
      \code{getTable(object)}: Retrieves the indicated table as a
      \code{data.frame}. Note that not all tables are output in
      parseable form, and that UCSC will truncate responses if they
      exceed certain limits (usually around 100,000 records). The safest
      (and most efficient) bet for large queries is to download the file
      via FTP and query it locally.
    }
    \item{}{
      \code{tableNames(object)}: Gets the names of the tables available
      for the provider, table and range specified by the query.
    }
  }
}

\section{Accessor methods}{
  In the code snippets below, \code{x}/\code{object} is a
  \code{UCSCTableQuery} object.

  \describe{
    \item{}{\code{genome(x)}, \code{genome(x) <- value}: Gets or sets
    the genome identifier (e.g. \dQuote{hg18}) of the object.
    }
    \item{}{\code{hubUrl(x)}, \code{hubUrl(x) <- value}: Gets or sets
    the TrackHub URI.
    }
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
      a \code{IntegerRangesList}.
    }
    \item{}{\code{names(x)}, \code{names(x) <- value}: Get or set the
      names of the features to retrieve. If \code{NULL}, this filter is
      disabled.
    }
    \item{}{\code{ucscSchema(x)}: Get
      the \code{\linkS4class{UCSCSchema}} object describing the selected table.}
    \item{}{\code{ucscTables(genome, track)}: Get the list of tables for the 
      specified track(e.g. \dQuote{Assembly}) and genome identifier (e.g. \dQuote{hg19}).
      Here \code{genome} and \code{track} must be a single non-NA string.
    }
  }
}

\author{ Michael Lawrence }

\examples{
\dontrun{
# query using `session` provider
session <- browserSession()
genome(session) <- "mm9"
## choose the phastCons30way table for a portion of mm9 chr1
query <- ucscTableQuery(session, table = "phastCons30way",
                        range = GRangesForUCSCGenome("mm9", "chr12",
                                             IRanges(57795963, 57815592)))
## list the table names
tableNames(query)
## retrieve the track data
track(query)  # a GRanges object
## get the multiz30waySummary track
tableName(query) <- "multiz30waySummary"
## get a data.frame summarizing the multiple alignment
getTable(query)

# query using `genome identifier` provider
query <- ucscTableQuery("hg18", table = "snp129",
                        names = c("rs10003974", "rs10087355", "rs10075230"))
ucscSchema(query)
getTable(query)

# query using `TrackHub URI` provider
query <- ucscTableQuery("https://ftp.ncbi.nlm.nih.gov/snp/population_frequency/TrackHub/20200227123210/",
                        genome = "hg19", table = "ALFA_GLB")
getTable(query)
# get the list of tables for 'Assembly' track and 'hg19' genome identifier
ucscTables("hg19", "Assembly")
}
}

\keyword{methods}
\keyword{classes}
