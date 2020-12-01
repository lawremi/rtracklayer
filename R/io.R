
setClass("RTLFile", representation(resource = "character_OR_connection"),
         contains = "VIRTUAL")

setMethod("initialize", "RTLFile", function(.Object, ...) {
    .Deprecated("BiocFile", msg = "This class is extending the deprecated RTLFile class from rtracklayer. Use BiocFile from BiocIO in place of RTLFile.")
    callNextMethod()
})

setClass("CompressedFile", contains = c("RTLFile", "VIRTUAL"))

setMethod("initialize", "CompressedFile", function(.Object, ...) {
    .Deprecated("CompressedFile", msg = "This class is extending the deprecated CompressedFile class from rtracklayer. Use CompressedFile from BiocIO in place of CompressedFile from rtracklayer.")
    callNextMethod()
})

setClass("RTLFileList",
         prototype = prototype(elementType = "RTLFile"),
         contains = "SimpleList")

RTLFileList <- function(files) {
    new("RTLFileList", listData = files)
}

setMethod("showAsCell", "RTLFileList", function(object) {
    showAsCell(vapply(object, path, character(1L)))
})

setMethod("fileFormat", "RTLFile", function(x)
    tolower(sub("File$", "", class(x))))

setMethod("show", "RTLFile", function(object) {
  r <- resource(object)
  if (!isSingleString(r))
    r <- summary(r)$description
  cat(class(object), "object\nresource:", r, "\n")
})

setMethod("as.character", "RTLFile", function(x) path(x))

#resource <- function(x) {
#    .Deprecated("resource", msg = "Use BiocIO::resource()")
#    BiocIO::resource(x)
#}

path <- function(object) {
    .Deprecated("path", msg = "Use BiocIO::path")
    BiocIO::path(object)
}

FileForFormat <- function(path, format = file_ext(path)) {
    .Deprecated("FileForFormat", msg = "Use BiocIO::FileForFormat")
    BiocIO::FileForFormat(path, format)
}

.ConnectionManager <- setRefClass("ConnectionManager",
                                  fields = c(connections = "list"))

manager <- function() .ConnectionManager()

connection <- function(manager, x, open = "") {
  connectionForResource(manager, resource(x), open = open)
}

resourceDescription <- function(x) {
  r <- resource(x)
  if (is(r, "connection"))
    r <- summary(r)$description
  r
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

## First checks for Windows drive letter.
## There are no known URI schemes that are only a single character.
isURL <- function(uri) {
    if (!isSingleString(uri))
        return(FALSE)
    windowsDriveLetter <- .Platform$OS.type == "windows" &&
        grepl("^[A-Za-z]:[/\\]", uri)
    grepl("^[A-Za-z]+:", uri) && !windowsDriveLetter
}

## Uses XML::parseURI, except custom check for whether it is a URL
.parseURI <- function(uri) {
  if (!isURL(uri)) {
    parsed <- parseURI("")
    parsed$path <- uri
  } else {
    parsed <- parseURI(uri)
    if (parsed$scheme == "file" && .Platform$OS.type == "windows") 
      parsed$path <- substring(parsed$path, 2) # trim '/' from '/C:/foo/bar.txt'
  }
  parsed
}

normURI <- function(x) {
  if (!isSingleString(x))
    stop("URI must be a single, non-NA string")
  uri <- .parseURI(x)
  if (uri$scheme == "") # /// (vs. //) needed for Windows
    x <- paste("/", file_path_as_absolute(x), sep = "")
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
  if (uriIsLocal(x)) {
    exists <- file.exists(uri$path)
  } else {
    txt <- getURL(x, header = TRUE)
    exists <- grepl("^HTTP/\\d+\\.\\d+ 200 OK", txt)
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

connectionForResource <- function(manager, x, open = "") {
  resource <- decompress(manager, x)
  if (is.character(resource)) {
    if (!nzchar(resource))
      stop("path cannot be an empty string")
    uri <- .parseURI(resource)
    if (uri$scheme != "")
      con <- url(resource)
    else con <- file(resource)
  } else con <- resource
  if (!isOpen(con) && nzchar(open)) {
      open(con, open)
      con <- manage(manager, con)
  }
  con
}

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

manage <- function(manager, con) {
    manager$connections <- unique(c(manager$connections, list(con)))
    attr(con, "manager") <- manager
    con
}

managed <- function(manager, con) {
    con %in% manager$connections
}

unmanage <- function(manager, con) {
    manager$connections <- setdiff(manager$connections, con)
    attr(con, "manager") <- NULL
    con
}

release <- function(manager, con) {
    if (managed(manager, con)) {
        unmanage(manager, con)
        close(con)
    }
    con
}
