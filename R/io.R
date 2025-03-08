manager <- BiocIO:::manager

connection <- BiocIO:::connection

resourceDescription <- BiocIO:::resourceDescription

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

## First checks for Windows drive letter.
## There are no known URI schemes that are only a single character.
isURL <- BiocIO:::isURL

.parseURI <- BiocIO:::.parseURI

checkURI <- function(x) {
  if (!isSingleString(x))
    stop("URI must be a single, non-NA string")
  x
}

createResource <- function(x, dir = FALSE, content = "") {
  uri <- .parseURI(x)
  if (uri$scheme == "file" || uri$scheme == "") {
    if (!file.exists(uri$path)) {
      if (dir)
        dir.create(uri$path, recursive = TRUE)
      else writeLines(content, uri$path)
    } else warning("Path '", uri$path, "' already exists")
  } else stop("Cannot create a resource that is not a local file")
}

uriExists <- function(x) {
  uri <- .parseURI(x)
  if (uriIsLocal(uri)) {
    exists <- file.exists(uri$path)
  } else {
    exists <- HEAD(x)$status_code == 200L
  }
  exists
}

uriIsLocal <- function(x) {
  x$scheme == "file" || x$scheme == ""
}

uriIsWritable <- function(x) {
  uri <- .parseURI(x)
  if (uriIsLocal(uri)) {
    !file.access(uri$path, 2) ||
    (!file.exists(uri$path) && uriIsWritable(dirname(uri$path)))
  } else FALSE
}

checkArgFormat <- function(con, format) {
  if (toupper(format) !=
      substring(toupper(sub("File$", "", class(con))), 1, nchar(format)))
    stop("Cannot treat a '", class(con), "' as format '", format, "'")
}

connectionForResource <- BiocIO:::connectionForResource

## BestFileFormat

setGeneric("bestFileFormat",
           function(x, dest, ...) standardGeneric("bestFileFormat"))

setMethod("bestFileFormat", c("GenomicRanges", "ANY"),
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

setMethod("bestFileFormat", c("IntegerRangesList", "ANY"), function(x, dest) {
  "bed" # just ranges...
})

## Connection management (similar to memory management)

manage <- BiocIO:::manage

managed <- BiocIO:::managed

unmanage <- BiocIO:::unmanage

release <- BiocIO:::release
