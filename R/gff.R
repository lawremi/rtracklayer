# input and output for the GFF format

setGeneric("export.gff",
           function(object, con, version = c("1", "2", "3"),
                    source = "rtracklayer", append = FALSE, ...)
           standardGeneric("export.gff"))

setMethod("export.gff", "ANY",
          function(object, con, version, source, append, ...)
          {
            cl <- class(object)
            object <- try(as(object, "RangedData"), silent = TRUE)
            if (class(object) == "try-error")
              stop("cannot export object of class '", cl, "'")
            export.gff(object, con=con, version=version, source=source,
                       append=append, ...)
          })

setMethod("export.gff", c("RangedData", "characterORconnection"),
          function(object, con, version = c("1", "2", "3"), source, append,
                   index = FALSE)
{
  version <- match.arg(version)
  
  if (!append) {
    cat("", file = con) # clear any existing file
    gffComment(con, "gff-version", version)
    sourceVersion <- try(package.version(source), TRUE)
    if (!inherits(sourceVersion, "try-error"))
      gffComment(con, "source-version", source, sourceVersion)
    gffComment(con, "date", format(Sys.time(), "%Y-%m-%d"))
  }

  if (index)
    object <- sortBySeqnameAndStart(object)

  seqname <- seqnames(object)
  if (is.null(object$ID))
    object$ID <- rownames(object)
  if (version == "3")
    seqname <- urlEncode(seqname, "a-zA-Z0-9.:^*$@!+_?|-")
  if (!is.null(object$source))
    source <- object$source
  if (version == "3")
    source <- urlEncode(source, "\t\n\r;=%&,", FALSE)
  feature <- object$type
  if (is.null(feature))
    feature <- "sequence"
  score <- score(object)
  if (is.null(score)) {
    if (version == "1")
      score <- 0
    else score <- NA
  } else {
    if (!("score" %in% colnames(object)))
      colnames(object)[1] <- "score" ## avoid outputting as attribute
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
      attrs <- seqname
  } else {
    builtin <- c("type", "strand", "score", "phase", "source")
    custom <- setdiff(colnames(object), builtin)
    if (length(custom)) {
      if (version == "3") tvsep <- "=" else tvsep <- " "
      attrs <- unlist(values(object))
      attrs <- as.data.frame(sapply(custom, function(name) {
        x <- attrs[[name]]
        x_flat <- if (is(x, "List")) unlist(x, use.names=FALSE) else x
        x_char <- as.character(x_flat)
        x_char <- sub(" *$", "", sub("^ *", "", as.character(x_char))) # trim
        if (version == "3")
          x_char <- urlEncode(x_char, "\t\n\r;=%&,", FALSE)
        if (is(x, "List")) {
          x_char[is.na(x_char)] <- "."
          x_char <- pasteCollapse(relist(x_char, x))
          x_char[elementLengths(x) == 0] <- NA
        }
        x_char[is.na(x_char)] <- "\r"
        paste(name, x_char, sep = tvsep)
      }, simplify = FALSE))
      attrs <- do.call(paste, c(attrs, sep = ";"))
      attrs <- gsub("[^;]*?\r(;|$)", "", attrs)
      ## attrs <- apply(attrs, 1, function(row) {
      ##   paste(colnames(attrs)[!is.na(row)], row[!is.na(row)], sep = tvsep,
      ##         collapse="; ")
      ## })
      attrs[nchar(attrs) == 0] <- NA
    }
  }
  
  scipen <- getOption("scipen")
  options(scipen = 100) # prevent use of scientific notation
  on.exit(options(scipen = scipen))
  
  if (!is.null(attrs)) { # write out the rows with attributes first
    write.table(cbind(table, attrs)[!is.na(attrs),], con, sep = "\t", na = ".",
                quote = FALSE, col.names = FALSE, row.names = FALSE,
                append = TRUE)
    table <- table[is.na(attrs),]
  }
  
  write.table(table, con, sep = "\t", na = ".", quote = FALSE,
              col.names = FALSE, row.names = FALSE, append = TRUE)
  if (index)
    indexTrack(object, con, "gff")
  invisible(NULL)
})

setGeneric("import.gff",
           function(con, version = c("1", "2", "3"), genome = NA,
                    asRangedData = TRUE, colnames = NULL, ...)
           standardGeneric("import.gff"))
           
setMethod("import.gff", "characterORconnection",
          function(con, version = c("1", "2", "3"), genome,
                   asRangedData = TRUE, colnames = NULL)
{
  versionMissing <- missing(version)
  version <- match.arg(version)
  lines <- readLines(con, warn = FALSE) # unfortunately, not a table
  lines <- lines[nzchar(lines)]
  
  # check our version
  versionLine <- lines[grep("^##gff-version", lines)]
  if (length(versionLine)) {
    specVersion <- sub("^##gff-version *", "", versionLine)
    if (!versionMissing && specVersion != version)
      warning("gff-version directive indicates version is ", specVersion,
              ", not ", version)
    else version <- specVersion 
  }

  ## strip comments
  notComments <- grep("^[^#]", lines)
  lines <- lines[notComments]
  
### TODO: handle ontologies (store in RangedData)

  ## construct table
  fields <- c("seqname", "source", "type", "start", "end", "score", "strand",
              "phase", "attributes")
  linesSplit <- strsplit(lines, "\t", fixed=TRUE)
  fieldCounts <- elementLengths(linesSplit)
  if (any(fieldCounts > length(fields)) ||
      any(fieldCounts < (length(fields) - 1)))
    stop("GFF files must have ", length(fields), " tab-separated columns")
  haveAttr <- fieldCounts == length(fields)
  haveAttrMat <- do.call(rbind, linesSplit[haveAttr])
  noAttrMat <- do.call(rbind, linesSplit[!haveAttr])
  if (!is.null(noAttrMat))
    noAttrMat <- cbind(noAttrMat, "")
  table <- rbind(noAttrMat, haveAttrMat)
  colnames(table) <- fields
  
  # handle missings
  table[table == "."] <- NA

  attrCol <- table[,"attributes"]
  if (version == "3") {
    table <- table[,setdiff(colnames(table), "attributes"),drop=FALSE]
    table[table[,"strand"] == "?","strand"] <- NA
    tableDec <- urlDecode(as.vector(table))
    table <- matrix(tableDec, ncol=ncol(table), dimnames=dimnames(table))
  }

  extraCols <- c("type", "source", "phase", "strand")
  if (!is.null(colnames))
    extraCols <- intersect(extraCols, colnames)
  xd <- as(table[,extraCols,drop=FALSE], "DataFrame")

  if (is.null(colnames) || length(setdiff(colnames, c(extraCols, "score")))) {
    if (version == "1") {
      if (is.null(colnames) || "group" %in% colnames)
        attrList <- list(group = attrCol)
      else attrList <- list()
    } else {
      attrSplit <- strsplit(attrCol, ";")
      lines <- rep(seq_along(attrSplit), elementLengths(attrSplit))
      attrs <- sub(" *$", "", sub("^ *", "", unlist(attrSplit)))
      if (version == "3") {
        attrs <- paste(attrs, "=", sep = "")
        tvSplit <- strsplit(attrs, "=", fixed=TRUE)
        if (any(elementLengths(tvSplit) != 2))
          stop("Some attributes do not conform to the 'tag=value' format")
        tvMat <- matrix(unlist(tvSplit), nrow = 2)
        tags <- urlDecode(tvMat[1,])
        vals <- urlDecode(tvMat[2,])
      } else { # split on first space (FIXME: not sensitive to quotes)
        tags <- sub(" .*", "", attrs) # strip surrounding quotes
        vals <- sub("^\"([^\"]*)\"$", "\\1", sub("^[^ ]* ", "", attrs))
      }
      if (!is.null(colnames)) {
        keep <- tags %in% colnames
        lines <- lines[keep]
        vals <- vals[keep]
        tags <- tags[keep]
      }
      lineValByTag <- split.data.frame(cbind(lines, vals), tags)
      attrList <- lapply(lineValByTag, function(tag) {
        vals <- tag[,"vals"]
        if (any(grepl(",", vals))) {
          vals <- CharacterList(strsplit(vals, ",", fixed=TRUE))
          coerced <- suppressWarnings(as(vals, "NumericList"))
          if (!any(any(is.na(coerced))))
            vals <- coerced
          vec <- as(rep(list(character()), nrow(table)), class(vals))
        } else {
          coerced <- suppressWarnings(as.numeric(vals))
          if (!any(is.na(coerced)))
            vals <- coerced
          vec <- rep(NA, nrow(table))
        }
        vec[as.integer(tag[,"lines"])] <- vals
        vec
      })
    }
    xd <- DataFrame(xd, attrList)
  }
  suppressWarnings(score <- as.numeric(table[,"score"]))
  if (!all(is.na(score)) && (is.null(colnames) || "score" %in% colnames))
    xd$score <- score

  end <- as.integer(table[,"end"])
  GenomicData(IRanges(as.integer(table[,"start"]), end),
              xd, chrom = table[,"seqname"], genome = genome,
              asRangedData = asRangedData)
})

setGeneric("export.gff1",
           function(object, con, ...) standardGeneric("export.gff1"))
setMethod("export.gff1", "ANY",
          function(object, con, ...) export.gff(object, con, "1", ...))

setGeneric("export.gff2",
           function(object, con, ...) standardGeneric("export.gff2"))
setMethod("export.gff2", "ANY",
          function(object, con, ...) export.gff(object, con, "2", ...))

setGeneric("export.gff3",
           function(object, con, ...) standardGeneric("export.gff3"))
setMethod("export.gff3", "ANY",
          function(object, con, ...) export.gff(object, con, "3", ...))


setGeneric("import.gff1",
           function(con, ...) standardGeneric("import.gff1"))
setMethod("import.gff1", "ANY",
          function(con, ...) import.gff(con, "1", ...))

setGeneric("import.gff2",
           function(con, ...) standardGeneric("import.gff2"))
setMethod("import.gff2", "ANY",
          function(con, ...) import.gff(con, "2", ...))

setGeneric("import.gff3",
           function(con, ...) standardGeneric("import.gff3"))
setMethod("import.gff3", "ANY",
          function(con, ...) import.gff(con, "3", ...))


## Conversion to GFF-like structure

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
            values(children)$Parent <- rep(parentIds, elementLengths(x))
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

# utilities

gffComment <- function(con, ...)
    cat("##", paste(...), "\n", sep = "", file = con, append = TRUE)
