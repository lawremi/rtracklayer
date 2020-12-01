
.ConnectionManager <- setRefClass("ConnectionManager",
                                  fields = c(connections = "list"))

manager <- function() .ConnectionManager()

resource <- function(x) x@resource

`resource<-` <- function(x, value) {
    x@resource <- value
    x
}

connection <- function(manager, x, open = "") {
  connectionForResource(manager, resource(x), open = open)
}

resourceDescription <- function(x) {
  r <- resource(x)
  if (is(r, "connection"))
    r <- summary(r)$description
  r
}

FileForFormat <- function(path, format = file_ext(path)) {
  if (!(isSingleString(path) || is(path, "connection")))
    stop("'path' must be a single string or a connection object")
  if (!isSingleString(format))
    stop("'format' must be a single string")
  if (format == "")
    stop("Cannot detect format (no extension found in file name)")
  fileClassName <- paste0(format, "File")
  signatureClasses <- function(fun, pos) {
    matrix(unlist(findMethods(fun)@signatures), 3)[pos,]
  }
  fileClassNames <- unique(c(signatureClasses(export, 2),
                             signatureClasses(import, 1)))
  fileClassNames <- fileClassNames[grepl("File$", fileClassNames)]
  fileSubClassNames <- unlist(lapply(fileClassNames, function(x) {
    names(getClassDef(x)@subclasses)
  }), use.names = FALSE)
  fileClassNames <- c(fileClassNames, fileSubClassNames) 
  fileClassIndex <- match(tolower(fileClassName),
                          tolower(fileClassNames))
  if (is.na(fileClassIndex))
    stop("Format '", format, "' unsupported")
  fileClassName <- fileClassNames[fileClassIndex]
  fileClass <- getClass(fileClassName)
  pkg <- packageSlot(fileClass)
  if (is.null(pkg) || identical(pkg, ".GlobalEnv"))
    ns <- topenv()
  else ns <- getNamespace(pkg[1])
  constructorName <- fileClassName
  if(!exists(constructorName, ns)) {
    parentClassNames <- names(getClass(constructorName)@contains)
    constructorName <- names(which(sapply(parentClassNames, exists, ns)))[1]
    if (is.na(constructorName))
      stop("No constructor found for ", fileClassName)
  }
  get(constructorName, ns)(path)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

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

setMethod("bestFileFormat", c("RangesList", "ANY"), function(x, dest) {
  "bed" # just ranges...
})

## Uses XML::parseURI, except first checks for Windows drive letter.
## There are no known URI schemes that are only a single character.
.parseURI <- function(uri) {
  windowsDriveLetter <- .Platform$OS.type == "windows" &&
      grepl("^[A-Za-z]:[/\\]", uri)
  hasScheme <- grepl("^[A-Za-z]+:", uri) && !windowsDriveLetter
  if (!hasScheme) {
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
    x <- paste("file:///", file_path_as_absolute(x), sep = "")
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

## Connection management (similar to memory management)

manage <- function(con) {
  if (!is.null(attr(con, "finalizerEnv")))
    return(con)
  env <- new.env()
  finalizer <- function(obj) {
    if (exists("con", parent.env(environment()), inherits=FALSE)) {
      close(con)
      rm(con, inherits = TRUE)
      TRUE
    } else FALSE
  }
  env$finalizer <- finalizer
  reg.finalizer(env, finalizer)
  attr(con, "finalizerEnv") <- env
  rm(env)
  con
}

unmanage <- function(con) {
  attr(con, "finalizerEnv") <- NULL
  con
}

release <- function(con) {
  env <- attr(con, "finalizerEnv")
  if (!is.null(env))
    env$finalizer()
  else FALSE
}

