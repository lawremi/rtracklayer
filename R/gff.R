### =========================================================================
### GFF (General Feature Format) support (all three versions, plus GTF)
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Classes
###

setClass("GFFFile", contains = "BiocFile")

## private
GFFFile <- function(resource, version = c("", "1", "2", "3")) {
  version <- match.arg(version)
  new(gffFileClass(version), resource = resource)
}

setClass("GFF1File", contains = "GFFFile")
GFF1File <- function(resource) {
  GFFFile(resource, "1")
}

setClass("GFF2File", contains = "GFFFile")
GFF2File <- function(resource) {
  GFFFile(resource, "2")
}

setClass("GFF3File", contains = "GFFFile")
GFF3File <- function(resource) {
  GFFFile(resource, "3")
}

setClass("GTFFile", contains = "GFF2File")
GTFFile <- function(resource) {
  new("GTFFile", GFF2File(resource))
}

setClass("GVFFile", contains = "GFF3File")
GVFFile <- function(resource) {
  new("GVFFile", GFF3File(resource))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export
###

setGeneric("export.gff",
           function(object, con, ...) standardGeneric("export.gff"))

setMethod("export.gff", "ANY",
          function(object, con, ...)
          {
            export(object, con, ...)
          })

setMethod("export", c("ANY", "GFFFile"),
          function(object, con, format, ...)
          {
            if (hasMethod("asGFF", class(object)))
              object <- asGFF(object)
            res <- try(as(object, "GRanges"), silent = TRUE)
            if (is(res, "try-error")) {
              res <- try(as(object, "GenomicRangedDataList"), silent = TRUE)
              if (is(res, "try-error"))
                stop("cannot export object of class '", class(object), "'")
            }
            object <- res
            if (!missing(format))
              checkArgFormat(con, format)
            export(object, con, ...)
          })

setMethod("export", c("CompressedGRangesList", "GFFFile"),
          function(object, con, format, ...)
          {
            object <- asGFF(object)
            callGeneric()
          }
          )

setMethod("export", c("GRangesList", "GTFFile"),
          function(object, con, format, ...) {
              stop("export of GRangesList to GTF is not yet supported")
              ## there is a start on asGTF() later in this file
              ## object <- asGTF(object)
              ## callGeneric()
          }
          )

setMethod("export", c("GenomicRanges", "GFFFile"),
          function(object, con, format, version = c("1", "2", "3"),
                   source = "rtracklayer", append = FALSE, index = FALSE)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            if (!missing(version) || !length(gffFileVersion(con)))
              con <- asGFFVersion(con, match.arg(version))
            version <- gffFileVersion(con)
            
            file <- con
            con <- resource(con)
            
            if (!append) {
              cat("", file = con) # clear any existing file
              gffComment(con, "gff-version", version)
              sourceVersion <- try(packageVersion(source), TRUE)
              if (!inherits(sourceVersion, "try-error"))
                gffComment(con, "source-version", source, sourceVersion)
              gffComment(con, "date", base::format(Sys.time(), "%Y-%m-%d"))
              genome <- singleGenome(genome(object))
              if (!is.na(genome))
                gffComment(con, "genome-build", paste(".", genome, sep = "\t"))
            }
            
            if (index)
              object <- sortBySeqnameAndStart(object)

            seqname <- seqnames(object)
            if (is.null(mcols(object)$ID))
              mcols(object)$ID <- names(object)
            if (version == "3")
              seqname <- urlEncode(seqname, "a-zA-Z0-9.:^*$@!+_?|-")
            if (!is.null(mcols(object)$source) && missing(source))
              source <- mcols(object)$source
            else source <- rep(source, length(object))
            if (version == "3")
              source <- urlEncode(source, "\t\n\r;=%&,", FALSE)
            feature <- mcols(object)$type
            if (is.null(feature))
              feature <- rep("sequence_feature", length(object))
            score <- score(object)
            if (is.null(score)) {
              score <- rep(NA_real_, length(object))
            } else {
              if (!("score" %in% colnames(mcols(object))))
                ## avoid outputting as attribute
                colnames(mcols(object))[1] <- "score" 
            }
            strand <- strand(object)
            if (is.null(strand))
                strand <- rep(strand(NA_character_), length(object))
            strand[strand == "*"] <- NA_integer_
            frame <- mcols(object)$phase
            if (is.null(frame)) {
              frame <- rep(NA_integer_, length(object))
              if ("CDS" %in% feature)
                warning(wmsg("The phase information is missing. ",
                             "The written file will contain CDS with ",
                             "no phase information."))
            } else {
              if (anyNA(frame[feature %in% "CDS"]))
                warning(wmsg("The phase information is missing for some CDS. ",
                             "The written file will contain some CDS with ",
                             "no phase information."))
            }
            
            table <- data.frame(seqname, source, feature, start(object),
                                end(object), score, strand, frame)

            attrs <- NULL
            if (version == "1") {
              attrs <- mcols(object)$group
              if (is.null(attrs))
                attrs <- as.vector(seqname)
            } else {
              builtin <- c("type", "score", "phase", "source")
              custom <- setdiff(colnames(mcols(object)), builtin)
              if (length(custom)) {
                if (version == "3") tvsep <- "=" else tvsep <- " "
                attrs <- mcols(object)
                attrs <- as.data.frame(sapply(custom, function(name) {
                  x <- attrs[[name]]
                  x_flat <- if (is(x, "List")) unlist(x, use.names=FALSE) else x
                  x_char <- as.character(x_flat)
                  x_char <- sub(" *$", "", sub("^ *", "", as.character(x_char)))
                  if (version == "3")
                    x_char <- urlEncode(x_char, "%\t\n\r;=&,", FALSE)
                  if (is(x, "List")) {
                    x_char[is.na(x_char)] <- "."
                    x_char <- pasteCollapse(relist(x_char, x))
                    x_char[elementNROWS(x) == 0] <- NA
                  }
                  ## FIXME: add option so these become "." instead of removing
                  x_char[is.na(x_char)] <- "\r"
                  if (!is.numeric(x_flat) && version != "3")
                    x_char <- paste0("\"", x_char, "\"")
                  paste(name, x_char, sep = tvsep)
                }, simplify = FALSE))
                if (version == "3") sep <- ";" else sep <- "; "
                attrs <- do.call(paste, c(attrs, sep = sep))
                attrs <- gsub("[^;]*?\r\"?(;|$)", "", attrs)
                attrs[nchar(attrs) == 0] <- NA
              }
            }
            
            scipen <- getOption("scipen")
            options(scipen = 100) # prevent use of scientific notation
            on.exit(options(scipen = scipen))
            
            if (!is.null(attrs)) { # write out the rows with attributes first
              write.table(cbind(table, attrs)[!is.na(attrs),], con, sep = "\t",
                          na = ".", quote = FALSE, col.names = FALSE,
                          row.names = FALSE, append = TRUE)
              table <- table[is.na(attrs),]
            }
            
            write.table(table, con, sep = "\t", na = ".", quote = FALSE,
                        col.names = FALSE, row.names = FALSE, append = TRUE)
            if (index)
              indexTrack(file)
            invisible(NULL)
          })

setMethod("export", c("SimpleGRangesList", "GFFFile"),
          .export_SimpleGRangesList_BiocFile)

setGeneric("export.gff1",
           function(object, con, ...) standardGeneric("export.gff1"))
setMethod("export.gff1", "ANY",
          function(object, con, ...) export(object, con, "gff1", ...))

setGeneric("export.gff2",
           function(object, con, ...) standardGeneric("export.gff2"))
setMethod("export.gff2", "ANY",
          function(object, con, ...) export(object, con, "gff2", ...))

setGeneric("export.gff3",
           function(object, con, ...) standardGeneric("export.gff3"))
setMethod("export.gff3", "ANY",
          function(object, con, ...) export(object, con, "gff3", ...))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import
###

setGeneric("import.gff", function(con, ...) standardGeneric("import.gff"))

setMethod("import.gff", "ANY",
          function(con, ...)
          {
            import(con, "gff", ...)
          })

setMethod("import", "GFFFile",
          function(con, format, text, version = c("", "1", "2", "3"),
                   genome = NA, colnames = NULL,
                   which = NULL, feature.type = NULL,
                   sequenceRegionsAsSeqinfo = FALSE)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            if (!missing(version))
              con <- asGFFVersion(con, match.arg(version))
            stopifnot(isTRUEorFALSE(sequenceRegionsAsSeqinfo))
            
            ## download the file first if it's remote
            if (is.character(resource(con))) {
                uri <- .parseURI(resource(con))
                if (uri$scheme %in% c("ftp", "http")) {
                    destfile <- tempfile()
                    download.file(resource(con), destfile)
                    con@resource <- destfile
                }
            }

            m <- manager()
            sniff_con <- connection(m, con, "r")
            on.exit(release(m, sniff_con))
            sniffed <- .sniffGFFVersion(sniff_con)
            version <- gffFileVersion(con)
            if (!length(version)) {
              if (is.null(sniffed))
                sniffed <- "1"
              con <- asGFFVersion(con, sniffed)
            }
            
            if (length(version) && !is.null(sniffed) &&
                !identical(sniffed, version))
              warning("gff-version directive indicates version is ", sniffed,
                      ", not ", version)

            if (is(genome, "Seqinfo") && length(genome) == 0L) {
                genome <- NA_character_
            }
               
            if (is.na(genome)) {
              genome <- genome(con)
              if (is.null(genome))
                  genome <- NA
            }
            
### FIXME: a queryForLines() function would be more efficient

            ## Temporarily disable use of Tabix Index.
            ## TODO: Restore use of Tabix Index!
            #con <- queryForResource(m, con, which)
            resource <- queryForResource(m, con)
            on.exit(release(m, resource), add=TRUE)
            ans <- readGFFAsGRanges(resource,
                                    version=version,
                                    colnames=colnames,
                                    filter=list(type=feature.type),
                                    genome=genome,
                                    sequenceRegionsAsSeqinfo=
                                        sequenceRegionsAsSeqinfo,
                                    speciesAsMetadata=TRUE)
            if (!attr(resource, "usedWhich") && !is.null(which))
                ans <- subsetByOverlaps(ans, which)
            ans
          })

setGeneric("import.gff1",
           function(con, ...) standardGeneric("import.gff1"))
setMethod("import.gff1", "ANY",
          function(con, ...) import(con, "gff1", ...))

setGeneric("import.gff2",
           function(con, ...) standardGeneric("import.gff2"))
setMethod("import.gff2", "ANY",
          function(con, ...) import(con, "gff2", ...))

setGeneric("import.gff3",
           function(con, ...) standardGeneric("import.gff3"))
setMethod("import.gff3", "ANY",
          function(con, ...) import(con, "gff3", ...))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DNAStringSet from fasta data
###



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setGeneric("asGFF", function(x, ...) standardGeneric("asGFF"))

setMethod("asGFF", "GRangesList",
          function(x, parentType = "mRNA", childType = "exon") {
            parent_range <- range(x)
            if (!all(elementNROWS(parent_range) == 1))
              stop("Elements in a group must be on same sequence and strand")
            parents <- unlist(parent_range, use.names = FALSE)
            children <- unlist(x, use.names = FALSE)
            makeId <- function(x, prefix) {
                paste(prefix, seq_len(length(x)), sep = "")
            }
            parentIds <- makeId(parents, parentType)
            values(parents)$type <- parentType
            values(parents)$ID <- parentIds
            values(parents)$Name <- names(x)
            values(children)$type <- childType
            values(children)$ID <- makeId(children, childType)
            values(children)$Name <- names(children)
            values(children)$Parent <- rep.int(parentIds, elementNROWS(x))
            allColumns <- union(colnames(values(parents)),
                                colnames(values(children)))
            values(children) <- rectifyDataFrame(values(children), allColumns)
            values(parents) <- rectifyDataFrame(values(parents), allColumns)
            c(parents, children)
          })

rectifyDataFrame <- function(x, allColumns) {
    x[setdiff(allColumns, colnames(x))] <- DataFrame(NA)
    x[allColumns]
}

### FIXME: We wrote this but never tested it, and it is not yet
### used. People should use GFF3 instead of this.
###
### KNOWN ISSUES:
### 1) The stop codon should not be included in the CDS (annoying)
### 2) pmapFromTranscripts() does not yet support our usage of it
### 3) Needs to move to GenomicFeatures

setGeneric("asGTF", function(x, ...) standardGeneric("asGTF"))

frame <- function(x) {
    cs <- cumsum(width(x))
    ucs <- unlist(cs, use.names=FALSE)
    ucs[end(PartitioningByEnd(x))] <- 0L
    ucs <- c(0L, head(ucs, -1L))
    ucs %% 3L
}

setMethod("asGTF", "GRangesList",
          function(x) {
              tx_ids <- names(x)
              if (is.null(tx_ids)) {
                  tx_ids <- seq_along(x)
              }
              processFeatures <- function(f) {
                  ans <- unlist(f, use.names=FALSE)
                  ans$frame <- frame(f)
                  if (is.null(ans$gene_id)) {
                      ans$gene_id <- ""
                  }
                  if (is.null(ans$transcript_id)) {
                      ans$transcript_id <- tx_ids[togroup(f)]
                  }
                  ans
              }
              start_codon_tx <-
                  GenomicFeatures::pmapFromTranscripts(IRanges(1L, 3L), x)
              start_codon <- processFeatures(start_codon_tx)
              mcols(start_codon)$type <- "start_codon"
              stop_ranges <- IRanges(end=sum(width(x)), width=3L)
              stop_codon_tx <-
                  GenomicFeatures::pmapFromTranscripts(stop_ranges, x)
              stop_codon <- processFeatures(stop_codon_tx)
              mcols(stop_codon)$type <- "stop_codon"
              codons <- c(start_codon, stop_codon)
              cds <- processFeatures(x)
              mcols(cds)$type <- "CDS"
              values(codons) <- rectifyDataFrame(values(codons), colnames(cds))
              c(codons, cds)
          })

## setMethod("asGTF", "TxDb",
##           function(x, by) {
##               cds <- cds(x, columns="tx_id")
##               cds <- cds[togroup(cds$tx_id)]
##               cds$transcript_id <- unlist(cds$tx_id)
##               cds$tx_id <- NULL
##               txGene <- transcriptsBy(x)
##               txToGene <- setNames(names(txGene)[togroup(txGene)],
##                                    unlist(txGene, use.names=FALSE)$tx_id)
##               cds$gene_id <- txToGene[cds$transcript_id]
##               processUTRs <- function(utr) {
##                   ans <- unlist(utr, use.names=FALSE)
##                   ans$transcript_id <- names(utr)[togroup(utr)]
##                   ans$gene_id <- txToGene[ans$transcript_id]
##                   ans
##               }
##               three_utr <- processUTRs(threeUTRsByTranscript(x))
##               three_utr$type <- "3UTR"
##               five_utr <- processUTRs(fiveUTRsByTranscript(x))
##               five_utr$type <- "5UTR"
##               c(asGTF(split(cds, cds$transcript_id)), five_utr, three_utr)
##           })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

scanGFFDirectives <- function(con, tag = NULL) {
  m <- manager()
  con <- connection(m, con, "r")
  on.exit(release(m, con))
  directives <- character()
  lines <- line <- readLines(con, n = 1)
  while(grepl("^#", line)) {
    if (grepl("^##", line)) {
        directives <- c(directives, line)
    }
    line <- readLines(con, n = 1)
    if (length(line) == 0L)
        break
    lines <- c(lines, line)
  }
  pushBack(lines, con)
  sub("^[^[:space:]]* ", "", grep(paste0("^##", tag), directives, value = TRUE))
}

gffGenomeBuild <- function(x) {
  genome_build <- scanGFFDirectives(x, "genome-build")
  unlist(strsplit(genome_build, "\t", fixed = TRUE))
}

setMethod("provider", "GFFFile", function(x) {
  gffGenomeBuild(x)[1]
})

setMethod("providerVersion", "GFFFile", function(x) {
  gffGenomeBuild(x)[2]
})
setMethod("genome", "GFFFile", function(x) providerVersion(x))

gffComment <- function(con, ...) 
  cat("##", paste(...), "\n", sep = "", file = con, append = TRUE)

.sniffGFFVersion <- function(con) {
  version <- NULL
  lines <- line <- readLines(con, n = 1)
  while(grepl("^#", line)) {
    if (grepl("^##gff-version", line)) {
      version <- sub("^##gff-version *", "", line)
      break
    }
    line <- readLines(con, n = 1)
    lines <- c(lines, line)
  }
  pushBack(lines, con)
  version
}

gffFileClass <- function(version) {
  paste("GFF", version, "File", sep = "")
}

gffFileVersion <- function(file) {
  versions <- c("1", "2", "3")
  unlist(Filter(function(v) is(file, gffFileClass(v)), versions))
}

asGFFVersion <- function(con, version) {
  if (!is(con, gffFileClass(version))) {
    if (class(con) != "GFFFile")
      warning("Treating a '", class(con), "' as GFF version '", version, "'")
    con <- GFFFile(resource(con), version)
  }
  con
}
