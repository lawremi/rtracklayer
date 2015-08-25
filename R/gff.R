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

### Return the 9 standard GFF columns as specified at:
###   http://www.sequenceontology.org/resources/gff3.html
GFFcolnames <- function() .Call(gff_colnames)

### 'df' must be a data-frame-like object (typically an ordinary data frame or
### a DataFrame object). 'decode_idx' must be a non-empty integer vector
### indicating which columns to decode. The columns to decode must be character
### vectors.
.urlDecodeCols <- function(df, decode_idx)
{
    decoded_cols <- lapply(setNames(decode_idx, colnames(df)[decode_idx]),
                           function(j)
                             urlDecode(df[[j]], na.strings=NA_character_))
    df[decode_idx] <- decoded_cols
    df
}

### 'df' must be a data-frame-like object (typically an ordinary data frame or
### a DataFrame object). 'split_idx' must be a non-empty integer vector
### indicating which columns to split. The columns to split must be character
### vectors. Split values are passed thru urlDecode() unless 'raw_data' is
### TRUE. Always returns a DataFrame.
.strsplitCols <- function(df, split_idx, raw_data)
{
    split_cols <- lapply(setNames(split_idx, colnames(df)[split_idx]),
        function(j) {
            col <- df[[j]]
            ## Probably the most efficient way to create an empty CharacterList
            ## of arbitrary length.
            split_col <- relist(character(0),
                                PartitioningByEnd(rep.int(0L, length(col))))
            not_na <- !is.na(col)
            tmp <- strsplit(col[not_na], ",", fixed=TRUE)
            split_col[not_na] <- CharacterList(tmp)
            if (raw_data)
                return(split_col)
            relist(urlDecode(unlist(split_col)), split_col)
        })
    ## Surprisingly sticking the CharacterList cols back into 'df' works
    ## even if 'df' is an ordinary data frame!
    df[split_idx] <- split_cols
    ans <- DataFrame(df, check.names=FALSE)
    ## "show" method for DataFrame is broken if some colnames are the empty
    ## string so we rename this column (in our case, we know there can only
    ## be one).
    m <- match("", colnames(ans))
    if (!is.na(m))
        colnames(ans)[m] <- "__empty_tag__"
    ans
}

### Does NOT work on a connection object.
readGFF <- function(filepath, columns=NULL, tags=NULL,
                    filter=NULL, raw_data=FALSE)
{
    if (!isSingleString(filepath))
        stop(wmsg("'filepath' must be a single string"))
    if (!isTRUEorFALSE(raw_data))
        stop(wmsg("'raw_data' must be TRUE or FALSE"))
    filexp <- XVector:::open_input_files(filepath)[[1L]]

    ## Check 'tags'.
    if (!is.null(tags)) {
        if (!is.character(tags))
            stop(wmsg("'tags' must be NULL or character vector"))
        if (any(is.na(tags)) || anyDuplicated(tags))
            stop(wmsg("'tags' cannot contain NAs or duplicates"))
    }

    ## Prepare 'colmap'.
    GFF_colnames <- GFFcolnames()
    stopifnot(GFF_colnames[[length(GFF_colnames)]] == "attributes")
    if (is.null(columns)) {
        colmap <- seq_along(GFF_colnames)
        ## We don't load the "attributes" column unless the user requested no
        ## tags (i.e. by setting 'tags' to character(0)).
        if (!(is.character(tags) && length(tags) == 0L))
            colmap[[length(GFF_colnames)]] <- NA_integer_
    } else if (is.character(columns)) {
        if (!all(columns %in% GFF_colnames)) {
            in1string <- paste0(GFF_colnames, collapse=", ")
            stop(wmsg("'columns' must contain valid GFF columns. ",
                      "Valid GFF columns are: ", in1string))
        }
        if (anyDuplicated(columns))
            stop(wmsg("'columns' cannot contain duplicates"))
        colmap <- match(GFF_colnames, columns)
    } else {
        stop(wmsg("'columns' must be NULL or character vector"))
    }

    ## Normalize 'filter'.
    if (!is.null(filter)) {
        if (!is.list(filter))
            stop("'filter' must be NULL or a named list")
        filter_names <- names(filter)
        if (is.null(filter_names))
            stop("'filter' must have names")
        valid_filter_names <- head(GFF_colnames, n=-1L)
        if (!all(filter_names %in% valid_filter_names)) {
            in1string <- paste0(valid_filter_names, collapse=", ")
            stop(wmsg("The names on 'filter' must be valid GFF columns ",
                      "(excluding \"attributes\"). ",
                      "Valid names on 'filter': ", in1string))
        }
        if (anyDuplicated(filter_names))
            stop(wmsg("names on 'filter' must be unique"))
        filter <- filter[valid_filter_names]
    }

    ## Return 'ans' as an ordinary data frame.
    ans <- .Call(gff_read, filexp, colmap, tags, filter, raw_data)
    ncol0 <- attr(ans, "ncol0")
    ntag <- attr(ans, "ntag")
    raw_data <- attr(ans, "raw_data")  # should be the same as user-supplied
    pragmas <- attr(ans, "pragmas")

    ## Post-process main GFF cols.
    #factor_colnames <- c("seqid", "source", "type", "strand")
    factor_colnames <- c("seqid", "source", "type")
    m <- match(factor_colnames, head(colnames(ans), n=ncol0))
    m <- m[!is.na(m)]
    factor_cols <- lapply(setNames(m, colnames(ans)[m]),
                          function(j)
                            factor(ans[[j]], levels=unique(ans[[j]])))
    ans[m] <- factor_cols

    ## Post-process tags.
    if (ntag != 0L) {
        multi_tags <- c("Parent", "Alias", "Note", "Dbxref", "Ontology_term")
        is_multi_tag <- sapply(seq_len(ntag) + ncol0,
                               function(j)
                                 colnames(ans)[[j]] %in% multi_tags ||
                                   any(grepl(",", ans[[j]], fixed=TRUE)))
        if (!raw_data) {
            decode_idx <- which(!is_multi_tag) + ncol0
            if (length(decode_idx) != 0L)
                ans <- .urlDecodeCols(ans, decode_idx)
        }
        split_idx <- which(is_multi_tag) + ncol0
        if (length(split_idx) != 0L) {
            ## Returns 'ans' as a DataFrame.
            ans <- .strsplitCols(ans, split_idx, raw_data)
        }
    }

    ## 'ans' could have lost its readGFF-specific attributes (e.g. if it was
    ## turned into a DataFrame), so we restore them and cross our fingers that
    ## they won't clash with the DataFrame slots the day the internals of
    ## DataFrame objects happen to change (very unlikely though).
    if (is.null(attr(ans, "ncol0")))
        attr(ans, "ncol0") <- ncol0
    if (is.null(attr(ans, "ntag")))
        attr(ans, "ntag") <- ntag
    if (is.null(attr(ans, "raw_data")))
        attr(ans, "raw_data") <- raw_data
    if (is.null(attr(ans, "pragmas")))
        attr(ans, "pragmas") <- pragmas
    ans
}

### sequence-region => Seqinfo
.parseSequenceRegionsAsSeqinfo <- function(lines) {
    sr <- grep("##sequence-region", lines, value=TRUE)
    srcon <- file()
    on.exit(close(srcon))
    writeLines(sr, srcon)
    srt <- read.table(srcon, comment.char="",
                      colClasses=list(NULL, "character", "integer",
                          "integer"))
    if (any(srt[[2L]] != 1L)) {
        warning("One or more ##sequence-region directives do not start at 1. ",
                "The assumptions made by 'sequenceRegionsAsSeqinfo=TRUE' ",
                "have been violated.")
    }
    Seqinfo(srt[[1L]], srt[[3L]])
}

.parseSpeciesAsMetadata <- function(lines) {
    species <- unique(grep("##species", lines, fixed=TRUE, value=TRUE))
    if (length(species) > 1L) {
        stop("multiple species definitions found")
    }
    metadata <- list()
    if (length(species) == 1L) {
        species <- sub("##species ", "", species, fixed=TRUE)
        if (isNCBISpeciesURL(species)) {
            ncbiError <- function(e) {
                warning("failed to retrieve organism information from NCBI")
            }
            metadata <- tryCatch(metadataFromNCBI(species), error = ncbiError)
        }
    }
    metadata
}

### Return an ordinary data frame with the 9 columns specified at
###   http://www.sequenceontology.org/resources/gff3.html
.read_gff <- function(con, isGFF3File=TRUE, feature.type=NULL,
                      sequenceRegionsAsSeqinfo=FALSE, speciesAsMetadata=FALSE)
{
    lines <- readLines(con, warn = FALSE) # unfortunately, not a table
    lines <- lines[nzchar(lines)]

    comments <- substr(lines, start=1L, stop=1L) == "#"
    comment_lines <- lines[comments]
    if (sequenceRegionsAsSeqinfo)
        ans_seqinfo <- .parseSequenceRegionsAsSeqinfo(comment_lines)
    if (speciesAsMetadata)
        ans_metadata <- .parseSpeciesAsMetadata(comment_lines)
    lines <- lines[!comments]
    rm(comments)

    ### TODO: handle ontologies (store in RangedData)

    ## strip FASTA sequence
    fastaHeaders <- which(substr(lines, start=1L, stop=1L) == ">")
    if (length(fastaHeaders))
        lines <- head(lines, fastaHeaders[1] - 1)

    ## construct table (character matrix)
    fields <- c("seqid", "source", "type", "start", "end", "score",
                "strand", "phase", "attributes")
    linesSplit <- strsplit(lines, "\t", fixed=TRUE)
    fieldCounts <- elementLengths(linesSplit)
    if (any(fieldCounts > length(fields)) ||
        any(fieldCounts < (length(fields) - 1)))
      stop("GFF files must have ", length(fields), " tab-separated columns")
    haveAttr <- fieldCounts == length(fields)
    data <- unlist(linesSplit[haveAttr], use.names=FALSE)
    if (is.null(data))
        data <- character(0)
    haveAttrMat <- matrix(data, ncol=length(fields), byrow=TRUE)
    data <- unlist(linesSplit[!haveAttr], use.names=FALSE)
    if (is.null(data))
        data <- character(0)
    rm(linesSplit)
    noAttrMat <- matrix(data, ncol=length(fields)-1L, byrow=TRUE)
    noAttrMat <- cbind(noAttrMat, rep.int("", nrow(noAttrMat)))
    table <- rbind(noAttrMat, haveAttrMat)
    rm(haveAttrMat, noAttrMat)
    colnames(table) <- fields

    if (!is.null(feature.type))
        table <- table[table[,"type"] %in% feature.type,,drop=FALSE]

    table[table == "."] <- NA_character_

    if (isGFF3File) {
      table[table[, "strand"] == "?", "strand"] <- NA_character_
      table[ , -9L] <- urlDecode(table[ , -9L], na.strings=NA_character_)
    }

    ## turn 'table' (matrix) into data frame
    ans_seqid <- table[ , "seqid"]
    ans_seqid <- factor(ans_seqid, levels=unique(ans_seqid))
    ans_source <- table[ , "source"]
    ans_source <- factor(ans_source, levels=unique(ans_source))
    ans_type <- table[ , "type"]
    ans_type <- factor(ans_type, levels=unique(ans_type))
    ans_start <- as.integer(table[ , "start"])
    ans_end <- as.integer(table[ , "end"])
    suppressWarnings(ans_score <- as.numeric(table[ , "score"]))
    ans_strand <- table[ , "strand"]
    ans_strand <- factor(ans_strand, levels=unique(ans_strand))
    ans_phase <- as.integer(table[ , "phase"])
    ans_attributes <- table[ , "attributes"]
    ans <- data.frame(seqid=ans_seqid,
                      source=ans_source,
                      type=ans_type,
                      start=ans_start,
                      end=ans_end,
                      score=ans_score,
                      strand=ans_strand,
                      phase=ans_phase,
                      attributes=ans_attributes,
                      stringsAsFactors=FALSE)

    if (sequenceRegionsAsSeqinfo)
        attr(ans, "seqinfo") <- ans_seqinfo
    if (speciesAsMetadata)
        attr(ans, "metadata") <- ans_metadata
    ans
}

## A fast replacement for .read_gff() that doesn't work on a connection object.
.read_gff2 <- function(con, isGFF3File=TRUE, feature.type=NULL,
                      sequenceRegionsAsSeqinfo=FALSE, speciesAsMetadata=FALSE)
{
    ans <- readGFF(con, filter=list(type=feature.type))
    pragmas <- attr(ans, "pragmas")
    if (sequenceRegionsAsSeqinfo)
        attr(ans, "seqinfo") <- .parseSequenceRegionsAsSeqinfo(pragmas)
    if (speciesAsMetadata)
        attr(ans, "metadata") <- .parseSpeciesAsMetadata(pragmas)
    if (isGFF3File)
        ans[ , "source"] <- urlDecode(ans[ , "source"],
                                      na.strings=NA_character_)
    ans
}

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
                   which = NULL, feature.type = NULL,
                   sequenceRegionsAsSeqinfo = FALSE)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            if (!missing(version))
              con <- asGFFVersion(con, match.arg(version))
            asRangedData <- normarg_asRangedData(asRangedData, "import")
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

            df <- .read_gff(con,
                            isGFF3File=is(file, "GFF3File"),
                            feature.type=feature.type,
                            sequenceRegionsAsSeqinfo=
                              sequenceRegionsAsSeqinfo && is(file, "GFF3File"),
                            speciesAsMetadata=TRUE)
            if (!is.null(attr(df, "seqinfo")))
                attr(con, "seqinfo") <- attr(df, "seqinfo")

            extraCols <- c("source", "type", "score", "strand", "phase")
            if (!is.null(colnames))
              extraCols <- intersect(extraCols, c(colnames, "strand"))
            xd <- as(df[ , extraCols, drop=FALSE], "DataFrame")

            if (is.null(colnames) || length(setdiff(colnames, extraCols))) {
              attrCol <- df[ , "attributes"]
              attrList <- .parse_attrCol(attrCol, file, colnames)
              xd <- DataFrame(xd, attrList)
            }
 
            GenomicData(IRanges(df[ , "start"], df[ , "end"]),
                        xd, chrom = df[ , "seqid"], genome = genome,
                        seqinfo = attr(con, "seqinfo"),
                        which = if (attr(con, "usedWhich")) NULL else which,
                        metadata = attr(df, "metadata"))
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
