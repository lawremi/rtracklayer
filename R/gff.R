### =========================================================================
### GFF (General Feature Format) support (all three versions, plus GTF)
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Classes
###

setClass("GFFFile", contains = "RTLFile")

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
                stop("cannot export object of class '", class(object), "': ",
                     res)
            }
            object <- res
            if (!missing(format))
              checkArgFormat(con, format)
            export(object, con, ...)
          })

setMethod("export", c("GRangesList", "GFFFile"),
          function(object, con, format, ...)
          {
            object <- asGFF(object)
            callGeneric()
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
              sourceVersion <- try(package.version(source), TRUE)
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
            if (is.null(object$ID))
              object$ID <- names(object)
            if (version == "3")
              seqname <- urlEncode(seqname, "a-zA-Z0-9.:^*$@!+_?|-")
            if (!is.null(object$source) && missing(source))
              source <- object$source
            if (version == "3")
              source <- urlEncode(source, "\t\n\r;=%&,", FALSE)
            feature <- object$type
            if (is.null(feature))
              feature <- "sequence_feature"
            score <- score(object)
            if (is.null(score)) {
              score <- NA
            } else {
              if (!("score" %in% colnames(mcols(object))))
                ## avoid outputting as attribute
                colnames(mcols(object))[1] <- "score" 
            }
            strand <- strand(object)
            if (is.null(strand))
              strand <- NA
            frame <- object$phase
            if (is.null(frame))
              frame <- NA
            
            table <- data.frame(seqname, source, feature, start(object),
                                end(object), score, strand, frame)

            attrs <- NULL
            if (version == "1") {
              attrs <- object$group
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
                    x_char[elementLengths(x) == 0] <- NA
                  }
                  ## FIXME: add option so these become "." instead of removing
                  x_char[is.na(x_char)] <- "\r"
                  if (!is.numeric(x_flat) && version != "3")
                    x_char <- paste0("\"", x_char, "\"")
                  paste(name, x_char, sep = tvsep)
                }, simplify = FALSE))
                attrs <- do.call(paste, c(attrs, sep = "; "))
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

setMethod("export", c("GenomicRangesList", "GFFFile"),
          .export_GenomicRangesList_RTLFile)

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

.parse_attrCol <- function(attrCol, file, colnames) {
  if (is(file, "GFF1File")) {
    if (is.null(colnames) || "group" %in% colnames)
      attrList <- list(group = factor(attrCol, levels=unique(attrCol)))
    else attrList <- list()
  } else {
    if (length(attrCol) == 0L)
        return(list())
    sentinel <- "\b"
    con <- file()
    on.exit(close(con))
    writeLines(paste0(sub("; *$", "", attrCol), sentinel), con)
    attrs <- scan(con, what=character(),
                  strip.white = TRUE, quiet = TRUE, sep = ";",
                  quote = if (is(file, "GFF3File")) "" else "\"")
    lines <- togroup(PartitioningByEnd(grep(sentinel, attrs, fixed = TRUE)))
    attrs <- sub(sentinel, "", attrs, fixed = TRUE)
    tag.value.sep <- if (is(file, "GFF3File")) "=" else " "
    sep.pos <- regexpr(tag.value.sep, attrs, fixed=TRUE)
    if (any(sep.pos == -1L))
      stop("Some attributes do not conform to 'tag", tag.value.sep,
           "value' format")
    tags <- substring(attrs, 1L, sep.pos - 1L)
    vals <- substring(attrs, sep.pos + 1L, nchar(attrs))
    if (!is.null(colnames)) {
      keep <- tags %in% colnames
      lines <- lines[keep]
      vals <- vals[keep]
      tags <- urlDecode(tags[keep])
    }
    tags <- factor(tags, levels=unique(tags))
    lineByTag <- split(lines, tags)
    valByTag <- split(vals, tags)

    ## FIXME: Parent, Alias, Note, DBxref,
    ## Ontology_term are allowed to have multiple
    ## values. We should probably always return them as a
    ## CharacterList.
    multiTags <- c("Parent", "Alias", "Note", "DBxref", "Ontology_term")
    attrList <- sapply(levels(tags), function(tagName) {
      vals <- valByTag[[tagName]]
      if (is(file, "GFF3File") &&
          (any(grepl(",", vals, fixed=TRUE)) || tagName %in% multiTags)) {
        vals <- CharacterList(strsplit(vals, ",", fixed=TRUE))
        vals <- relist(urlDecode(unlist(vals)), vals)
        coerced <- suppressWarnings(as(vals, "NumericList"))
        if (!any(any(is.na(coerced))))
          vals <- coerced
        vec <- as(rep.int(list(character()), length(attrCol)), class(vals))
      } else {
        coerced <- suppressWarnings(as.numeric(vals))
        if (!any(is.na(coerced)))
          vals <- coerced
        if (is(file, "GFF3File"))
          vals <- urlDecode(vals)
        vec <- rep.int(NA, length(attrCol))
      }
      vec[lineByTag[[tagName]]] <- vals
      vec
    }, simplify = FALSE)
  }
  attrList
}

setMethod("import", "GFFFile",
          function(con, format, text, version = c("", "1", "2", "3"),
                   genome = NA, asRangedData = FALSE, colnames = NULL,
                   which = NULL, feature.type = NULL)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            if (!missing(version))
              con <- asGFFVersion(con, match.arg(version))
            asRangedData <- normarg_asRangedData(asRangedData, "import")

            sniffed <- sniffGFFVersion(resource(con))
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

            if (is.na(genome))
              genome <- genome(con)
            
### FIXME: a queryForLines() function would be more efficient
            file <- con
            con <- queryForResource(con, which)
            lines <- readLines(con, warn = FALSE) # unfortunately, not a table
            lines <- lines[nzchar(lines)]
            
            ## strip comments
            notComments <- which(substr(lines, start=1L, stop=1L) != "#")
            lines <- lines[notComments]
            
### TODO: handle ontologies (store in RangedData)

            ## strip FASTA sequence
            fastaHeaders <- which(substr(lines, start=1L, stop=1L) == ">")
            if (length(fastaHeaders))
              lines <- head(lines, fastaHeaders[1] - 1)

            ## construct table
            fields <- c("seqname", "source", "type", "start", "end", "score",
                        "strand", "phase", "attributes")
            linesSplit <- strsplit(lines, "\t", fixed=TRUE)
            fieldCounts <- elementLengths(linesSplit)
            if (any(fieldCounts > length(fields)) ||
                any(fieldCounts < (length(fields) - 1)))
              stop("GFF files must have ", length(fields),
                   " tab-separated columns")
            haveAttr <- fieldCounts == length(fields)
            data <- unlist(linesSplit[haveAttr], use.names=FALSE)
            if (is.null(data))
                data <- character(0)
            haveAttrMat <- matrix(data, ncol=length(fields), byrow=TRUE)
            data <- unlist(linesSplit[!haveAttr], use.names=FALSE)
            if (is.null(data))
                data <- character(0)
            noAttrMat <- matrix(data, ncol=length(fields)-1L, byrow=TRUE)
            noAttrMat <- cbind(noAttrMat, rep.int("", nrow(noAttrMat)))
            table <- rbind(noAttrMat, haveAttrMat)
            colnames(table) <- fields

            if (!is.null(feature.type))
                table <- table[table[,"type"] %in% feature.type,,drop=FALSE]

            ## handle missings
            table[table == "."] <- NA_character_
            
            attrCol <- table[,"attributes"]
            if (is(file, "GFF3File")) {
              table <- table[,setdiff(colnames(table), "attributes"),drop=FALSE]
              table[table[,"strand"] == "?","strand"] <- NA_character_
              is_not_NA <- !is.na(table)
              table[is_not_NA] <- urlDecode(table[is_not_NA])
            }
            
            extraCols <- c("source", "type", "score", "strand", "phase")
            if (!is.null(colnames))
              extraCols <- intersect(extraCols, colnames)
            xd <- as(table[,extraCols,drop=FALSE], "DataFrame")

            if (!is.null(xd$phase))
              xd$phase <- as.integer(as.character(xd$phase))
            if (!is.null(xd$score))
              suppressWarnings(xd$score <- as.numeric(as.character(xd$score)))

            if (is.null(colnames) || length(setdiff(colnames, extraCols))) {
              attrList <- .parse_attrCol(attrCol, file, colnames)
              xd <- DataFrame(xd, attrList)
            }
            
            end <- as.integer(table[,"end"])
            GenomicData(IRanges(as.integer(table[,"start"]), end),
                        xd, chrom = table[,"seqname"], genome = genome,
                        seqinfo = attr(con, "seqinfo"),
                        asRangedData = asRangedData,
                        which = if (attr(con, "usedWhich")) NULL else which)
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
            if (!all(elementLengths(parent_range) == 1))
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
            values(children)$Parent <- rep.int(parentIds, elementLengths(x))
            allColumns <- union(colnames(values(parents)),
                                colnames(values(children)))
            rectifyDataFrame <- function(x) {
              x[setdiff(allColumns, colnames(x))] <- DataFrame(NA)
              x[allColumns]
            }
            values(children) <- rectifyDataFrame(values(children))
            values(parents) <- rectifyDataFrame(values(parents))
            c(parents, children)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

scanGFFDirectives <- function(con, tag = NULL) {
  con <- connection(con, "r")
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

sniffGFFVersion <- function(con) {
  con <- connectionForResource(con, "r")
  version <- NULL
  lines <- line <- readLines(con, n = 1)
  while(grepl("^#", line)) {
    if (grepl("^##gff-version", line)) {
      version <- sub("^##gff-version ", "", line)
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
