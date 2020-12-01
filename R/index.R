### =========================================================================
### Separate index support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### queryForResource
###

## There is some question as to whether the import routines should
## recurse through import.tabix, or perform the query
## internally. Currently, we are taking the latter route, because it
## is consistent with (gzip) decoding and it allows the parsers more
## flexibility.

setGeneric("queryForResource",
           function(manager, x, which = NULL, ...)
               standardGeneric("queryForResource"),
           signature="x")

## Attaches 'usedWhich' attribute, an optimization hint indicating
## that subsetting by 'which' has been performed and is no longer
## necessary. Probably premature.

setMethod("queryForResource", "BiocFile", function(manager, x, which = NULL, ...)
{
  r <- resource(x)
  ans <- structure(r, usedWhich = FALSE)
  if (!is.null(which) && is.character(r)) {
    x_tbi <- paste(r, "tbi", sep = ".")
    if (file.exists(x_tbi))
      ans <- queryForResource(manager, TabixFile(r), which = which, ...)
  }
  ans
})

setMethod("queryForResource", "TabixFile",
          function(manager, x, which, header = TRUE, ...) {
            tabixHeader <- headerTabix(x)
            si <- Seqinfo(tabixHeader$seqnames)
            if (is.null(which)) {
              buffer <- connectionForResource(manager, path(x), "r")
              if (!header)
                readLines(buffer, tabixHeader$skip)
            } else {
              buffer <- manage(manager, file())
              if (header) {
                skippedLines <- readLines(path(x), tabixHeader$skip)
                writeLines(skippedLines, buffer)
              }
              lines <- unlist(scanTabix(x, param = which), use.names = FALSE)
              writeLines(lines, buffer)
              si <- merge(si, seqinfo(which))
            }
            structure(buffer, usedWhich = TRUE, seqinfo = si)
          })

queryForConnection <- function(manager, x, which = NULL, ...) {
  resource <- queryForResource(manager, x, which = which, ...)
  con <- connectionForResource(manager, resource, open = "r")
  structure(con, usedWhich = attr(resource, "usedWhich"))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### indexTrack: attempt to build an index for a track
###
### This is an obscene hack
###

indexTrack <- function(con, ...) {
  indexed <- NULL
  formats <- eval(formals(indexTabix)$format)
  uri <- path(con)
  parsed_uri <- .parseURI(uri)
  if (!uriIsLocal(parsed_uri))
    stop("'con' must be a path to a local file")
  original_path <- parsed_uri$path
  path <- bgzip(original_path, overwrite = TRUE)
  format <- Find(function(f) {
    is(con, paste(toupper(f), "File", sep = ""))
  }, formats)
  if (!is.null(format)) {
    indexTabix(path, format, ...)
  } else {
    indexTabix(path, ...)
  }
  indexed <- TabixFile(path)
  unlink(original_path)
  invisible(indexed)
}
