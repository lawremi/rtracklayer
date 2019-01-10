### =========================================================================
### Compression
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### General
###

setClass("CompressedFile", contains = c("RTLFile", "VIRTUAL"))

setGeneric("decompress",
           function(con, ...) standardGeneric("decompress"))

setMethod("decompress", "ANY", function(manager, con, ...) con,
          signature="con")

setMethod("decompress", "CompressedFile", function(manager, con, ...) {
  resource <- resource(con)
  if (is.character(resource))
    manage(manager, gzfile(resource)) # handles gzip, bzip2 and xz
  else stop("Cannot decompress connection")
})

setMethod("decompress", "character",
          function(manager, con, ...) {
            file <- try(FileForFormat(con), silent = TRUE)
            if (!is(file, "try-error")) {
              decompressed <- decompress(manager, file)
              if (!identical(file, decompressed))
                con <- decompressed
            }
            con
          })

## should only happen internally (user would not give compression as format)
setMethod("import", c("CompressedFile", "missing"),
          function(con, format, text, ...)
          {
            desc <- resourceDescription(con)
            con <- FileForFormat(resource(con),
                                 file_ext(file_path_sans_ext(desc)))
            import(con, ...)
          })

## 'compress' is a simple alias for 'decompress', since connections are two-way
compress <- decompress

## should only happen internally (user would not give compression as format)
setMethod("export", c("ANY", "CompressedFile", "missing"),
          function(object, con, format, ...)
          {
            desc <- resourceDescription(con)
            con <- FileForFormat(resource(con),
                                 file_ext(file_path_sans_ext(desc)))
            export(object, con, ...)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GZip
###

setClass("GZFile", contains = "CompressedFile")

GZFile <- function(resource) {
  new("GZFile", resource = resource)
}

setMethod("decompress", "GZFile", function(manager, con) {
  ungzip(manager, resource(con))
})

setGeneric("ungzip", function(manager, x, ...) standardGeneric("ungzip"),
           signature="x")

setMethod("ungzip", "character", function(manager, x) {
  uri <- .parseURI(x)
  if (uri$scheme != "" && uri$scheme != "file")
    con <- gzcon(url(x, open="rb"), text=TRUE)
  else con <- gzfile(uri$path)
  manage(manager, con)
})

setMethod("ungzip", "connection", function(x) {
  gzcon(x, text=TRUE)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### BGZip (Block-wise gzip files; extension used by Rsamtools::bgzip())
###

setClass("BGZFile", contains = "GZFile")

BGZFile <- function(resource) {
  new("BGZFile", resource = resource)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### BZip2
###

setClass("BZ2File", contains = "CompressedFile")

BZ2File <- function(resource) {
  new("BZ2File", resource = resource)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### XZ
###

setClass("XZFile", contains = "CompressedFile")

XZFile <- function(resource) {
  new("XZFile", resource = resource)
}
