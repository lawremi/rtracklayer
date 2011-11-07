### =========================================================================
### Import/export support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Classes files and connections
###

## Need to dispatch on S3 "connection"

##setOldClass("connection")

.connectionClasses <- c("file", "url", "gzfile", "bzfile", "unz", "pipe",
                        "fifo", "sockconn", "terminal", "textConnection",
                        "gzcon")
apply(cbind(.connectionClasses, "connection"), 1, setOldClass,
      where = environment())
setClassUnion("characterORconnection", c("character", "connection"))

setClass("RTLFile", representation(path = "character"), contains = "VIRTUAL")

setMethod("show", "RTLFile", function(object) {
  cat(class(object), "object\npath:", object@path, "\n")
})

setGeneric("path",
           function(object, ...) standardGeneric("path"))

setMethod("path", "RTLFile", function(object) object@path)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export
###

setGeneric("export",
           function(object, con, format, ...) standardGeneric("export"))

.exportForFormat <- function(format) {
  fun <- try(match.fun(paste("export", format, sep=".")), TRUE)
  if (is.character(fun))
    stop("No export function for '", format, "' found")
  fun
}

setMethod("export", c(con = "connection", format = "character"),
          function(object, con, format, ...)
          {
            if (!isOpen(con)) {
              open(con, "w")
              on.exit(close(con))
            }
            fun <- .exportForFormat(format)
            fun(object, con, ...)
          })
setMethod("export", c(con = "missing", format = "character"),
          function(object, con, format, ...)
          {
            con <- file()
            on.exit(close(con))
            export(object, con, format, ...)
            text <- readLines(con, warn = FALSE)
            text
          })
setMethod("export", c(con = "character", format = "missing"),
          function(object, con, format, ...)
          {
            ext <- file_ext(con)
            export(object, con, ext, ...)
          })
setMethod("export", c(con = "character", format = "character"),
          function(object, con, format, ...)
          {
            fun <- .exportForFormat(format)
            fun(object, con, ...)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import
###

.importForFormat <- function(format) {
  fun <- try(match.fun(paste("import", format, sep=".")), TRUE)
  if (is.character(fun))
    stop("No import function for '", format, "' found")
  fun
}

setGeneric("import",
           function(con, format, text, ...) standardGeneric("import"))

setMethod("import", c("connection", "character"),
          function(con, format, text, ...)
          {
            fun <- .importForFormat(format)
            if (!isOpen(con)) {
              open(con, "r")
              on.exit(close(con))
            }
            fun(con, ...)
          })
setMethod("import", c("character", "missing"),
          function(con, format, text, ...)
          {
            ext <- file_ext(con)
            import(con, ext, ...)
          })
setMethod("import", c("character", "character"),
          function(con, format, text, ...)
          {
            if (file_ext(con) == "gz") {
              if (format == "gz") # should only happen if user did not specify
                format <- file_ext(sub("\\.gz$", "", con))
              uri <- parseURI(con)
              if (uri$scheme != "")
                con <- gzcon(url(con))
              else con <- gzfile(con)
              import(con, format, ...)
            } else {
              fun <- .importForFormat(format)
              fun(con, ...)
            }
          })
setMethod("import", c(con = "missing", text = "character"),
          function(con, format, text, ...)
          {
            con <- file()
            writeLines(text, con)
            obj <- import(con, format, ...)
            close(con)
            obj
          })


## gzip handling
## setGeneric("import.gz",
##            function(con, ...) standardGeneric("import.gz"))

## setMethod("import.gz", "character", function(con, ...) {
##   import(gzfile(con), file_ext(sub("\\.gz$", "", con)), ...)
## })

## setMethod("import.gz", "connection", function(con, ...) {
##   import(gzcon(con), ...)
## })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

setGeneric("bestFileFormat",
           function(x, dest, ...) standardGeneric("bestFileFormat"))

setClassUnion("RangedDataORGenomicRanges", c("RangedData", "GenomicRanges"))

setMethod("bestFileFormat", c("RangedDataORGenomicRanges", "ANY"),
          function(x, dest) {
            ## have numbers on a single strand, use BigWig
            if (is.numeric(score(x)) && length(unique(strand(x))) == 1L)
              "bw"
            else "bed"
          })

setMethod("bestFileFormat", c("GRangesList", "ANY"), function(x, dest) {
  "bed" # need hierarchical structure
})

setMethod("bestFileFormat", c("RleList", "ANY"), function(x, dest) {
  "bw" # e.g., coverage
})

setMethod("bestFileFormat", c("RangesList", "ANY"), function(x, dest) {
  "bed" # just ranges...
})

file_ext <- function(con) gsub(".*\\.([^.]*)$", "\\1", con)

normURI <- function(x) {
  if (!isSingleString(x))
    stop("URI must be a single, non-NA string")
  uri <- parseURI(x)
  if (uri$scheme == "")
    x <- paste("file://", file_path_as_absolute(x), sep = "")
  x
}

createResource <- function(x, dir = FALSE, content = "") {
  uri <- parseURI(x)
  if (uri$scheme == "file" || uri$scheme == "") {
    if (!file.exists(uri$path)) {
      if (dir)
        dir.create(uri$path, recursive = TRUE)
      else writeLines(content, uri$path)
    } else warning("Path '", uri$path, "' already exists")
  } else stop("Cannot create a resource that is not a local file")
}

uriExists <- function(x) {
  uri <- parseURI(x)
  if (uriIsLocal(uri)) {
    exists <- file.exists(uri$path)
  } else {
    txt <- getURL(x, header = TRUE)
    exists <- grepl("^HTTP/\\d+\\.\\d+ 200 OK", txt)
  }
  exists
}

uriIsLocal <- function(x) {
  x$scheme == "file" || nchar(x$scheme) < 2 # allow for volumes on Windows
}

uriIsWritable <- function(x) {
  uri <- parseURI(x)
  if (uriIsLocal(uri)) {
    !file.access(uri$path, 2) ||
    (!file.exists(uri$path) && uriIsWritable(dirname(uri$path)))
  } else FALSE
}
