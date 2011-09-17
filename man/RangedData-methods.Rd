\name{RangedData-methods}
\docType{methods}
\alias{RangedData-methods}

% constructor
\alias{GenomicData}

% accessors
\alias{chrom}
\alias{chrom,RangedData-method}
\alias{chrom,GRanges-method}
\alias{chrom<-}
\alias{chrom<-,RangedData-method}
\alias{chrom<-,GRanges-method}
\alias{genome,RangedData-method}
%\alias{genome,GRanges-method}
\alias{genome<-,RangedData-method}
%\alias{genome<-,GRanges-method}

\title{Data on a Genome} \description{The \code{rtracklayer} package adds
convenience methods on top of \code{RangedData} and \code{GRanges} to
manipulate data on genomic ranges. For \code{RangedData} the spaces are
now called chromosomes (but could still refer to some other type of sequence).
Similarly the universe refers to the genome.}

\section{Accessors}{
  In the code snippets below,
  \code{x} is a \code{RangedData} or \code{GRanges} object.

  \describe{
    \item{}{
      \code{chrom(x), chrom(x) <- value}: Gets or sets the chromosome
      names for \code{x}. The length of \code{value} should equal the
      length of \code{x}. This is an alias for
      \code{\link[IRanges:RangedData-class]{names}(x)}.
    }
    \item{}{
      \code{genome(x)}, \code{genome(x) <- value}: Gets or sets the
      genome. If \code{x} is a \code{RangedData}, this is a single string
      or \code{NULL} and \code{genome}/\code{genome<-} are just simple
      wrappers around
      \code{\link[IRanges:RangedData-class]{universe}} and
      \code{\link[IRanges:RangedData-class]{universe<-}}, respectively.
      The \code{genome} and \code{genome<-} methods for \code{GRanges}
      objects are defined and described in the GenomicRanges package.
    }
  }
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenomicData(ranges, ..., strand = NULL, chrom = NULL,
        genome = NULL, asRangedData = TRUE)}: If \code{asRangedData} is
      \code{TRUE}, constructs a \code{RangedData} instance with the given
      \code{ranges} and variables in \code{...} (see the
      \code{\link[IRanges:RangedData-class]{RangedData}} constructor).
      If \code{asRangedData} is \code{FALSE}, constructs a \code{GRanges}
      instance with the given \code{ranges} and variables in \code{...}.

      If non-\code{NULL}, the \code{strand} argument specifies the strand
      of each range. It should be a character vector or factor of length
      equal to that of \code{ranges}. All values should be either \code{-},
      \code{+}, \code{*} or \code{NA}. (The \code{NA} code for \code{strand}
      is only acceptable when \code{asRangedData} is \code{TRUE}.) To get the
      levels for \code{strand}, call \code{levels(strand())}.

      \code{chrom} argument is analogous to \code{space} in the
      \code{RangedData} and \code{seqnames} in \code{GRanges} constructors.

      The \code{genome} argument should be a scalar string and is treated
      as the \code{RangedData} universe. See the examples.
      
      If \code{ranges} is not a \code{Ranges} object, this function
      calls \code{as(ranges, "RangedData")} and returns the result if
      successful. As a special case, the \dQuote{chrom} column in a
      \code{data.frame}-like object is renamed to \dQuote{space}, for
      convenience. Thus, one could pass a \code{data.frame} with columns
      \dQuote{start}, \dQuote{end} and, optionally, \dQuote{chrom}.
    }
  }
}

\author{ Michael Lawrence and Patrick Aboyoun }

\examples{
  range1 <- IRanges(start=c(1,2,3), end=c(5,2,8))

  ## just ranges ##
  ## RangedData instance
  rd <- GenomicData(range1)
  ## GRanges instance
  gr <- GenomicData(range1, asRangedData = FALSE)

  ## with a genome (universe) ##
  ## RangedData instance
  rd <- GenomicData(range1, genome = "hg18")
  genome(rd) ## "hg18"
  ## GRanges instance
  gr <- GenomicData(range1, genome = "hg18", asRangedData = FALSE)
  genome(gr) ## "hg18"

  ## with some data ##
  filter <- c(1L, 0L, 1L)
  score <- c(10L, 2L, NA)
  strand <- factor(c("+", NA, "-"), levels = levels(strand()))
  ## RangedData instance
  rd <- GenomicData(range1, score, genome = "hg18")
  rd[["score"]]
  strand(rd) ## all NA
  rd <- GenomicData(range1, score, filt = filter, strand = strand)
  rd[["filt"]]
  strand(rd) ## equal to 'strand'
  ## GRanges instance
  gr <- GenomicData(range1, score, genome = "hg18", asRangedData = FALSE)
  values(gr)[["score"]]
  strand(gr) ## all '*'
  gr <- GenomicData(range1, score, filt = filter, strand = strand,
                    asRangedData = FALSE)
  values(gr)[["filt"]]
  strand(gr) ## equal to 'strand'

  ## multiple chromosomes ##
  range2 <- IRanges(start=c(15,45,20,1), end=c(15,100,80,5))
  ranges <- c(range1, range2)
  score <- c(score, c(0L, 3L, NA, 22L)) 
  chrom <- paste("chr", rep(c(1,2), c(length(range1), length(range2))), sep="")
  ## RangedData instance
  rd <- GenomicData(ranges, score, chrom = chrom, genome = "hg18")
  chrom(rd) # equal to 'chrom'
  rd[["score"]] # unlists over the chromosomes
  score(rd)
  rd[1][["score"]] # equal to score[1:3]
  ## GRanges instance
  gr <- GenomicData(ranges, score, chrom = chrom, genome = "hg18",
                    asRangedData = FALSE)
  chrom(gr) # equal to 'chrom'
  values(gr)[["score"]]
  values(gr[chrom(gr) == "chr1"])[["score"]]

  ## coercion from data.frame ##
  df <- as.data.frame(rd)
  GenomicData(df)
  GenomicData(df, asRangedData = FALSE)
}
\keyword{classes}
\keyword{methods}
