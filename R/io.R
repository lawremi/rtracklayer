### =========================================================================
### Import/export support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Classes files and connections
###

### RTLFile is a base class for high-level file abstractions, where
### subclasses are associated with a particular file format/type. It
### wraps a low-level representation of a file, currently either a
### path/URL or connection.

setClass("RTLFile", representation(resource = "characterORconnection"),
         contains = "VIRTUAL")

resource <- function(x) x@resource

connection <- function(x, open = "") {
  connectionForResource(resource(x), open = open)
}

resourceDescription <- function(x) {
  r <- resource(x)
  if (is(r, "connection"))
    r <- summary(r)$description
  r
}

fileFormat <- function(x) {
  tolower(sub("File$", "", class(x)))
}

setMethod("path", "RTLFile", function(object) {
  r <- resource(object)
  if (!is.character(r))
    stop("Connection resource requested as a path")
  r
})

setMethod("show", "RTLFile", function(object) {
  r <- resource(object)
  if (!isSingleString(r))
    r <- summary(r)$description
  cat(class(object), "object\nresource:", r, "\n")
})

FileForFormat <- function(path, format = file_ext(path)) {
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
  if (is.na(fileClassIndex)) {
    stop("Format '", format, "' unsupported")
  }
  fileClassName <- fileClassNames[fileClassIndex]
  fileClass <- getClass(fileClassName)
  pkg <- packageSlot(fileClass)
  if (is.null(pkg))
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
### Export
###

setGeneric("export",
           function(object, con, format, ...) standardGeneric("export"))

setMethod("export", c(con = "connection", format = "character"),
          function(object, con, format, ...)
          {
            export(object, FileForFormat(con, format), ...)
          })

setMethod("export", c(con = "connection", format = "missing"),
          function(object, con, format, ...)
          {
            format <- file_ext(summary(con)$description)
            export(object, con, format, ...)
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
            export(object, FileForFormat(con), ...)
          })
setMethod("export", c(con = "character", format = "character"),
          function(object, con, format, ...)
          {
            export(object, FileForFormat(con, format), ...)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import
###

setGeneric("import",
           function(con, format, text, ...) standardGeneric("import"))

setMethod("import", c("connection", "character"),
          function(con, format, text, ...)
          {
            import(FileForFormat(con, format), ...)
          })
setMethod("import", c("connection", "missing"),
          function(con, format, text, ...)
          {
            format <- file_ext(summary(con)$description)
            import(con, format, ...)
          })
setMethod("import", c("character", "missing"),
          function(con, format, text, ...)
          {
            import(FileForFormat(con), ...)
          })
setMethod("import", c("character", "character"),
          function(con, format, text, ...)
          {
            import(FileForFormat(con, format), ...)
          })
setMethod("import", c(con = "missing", text = "character"),
          function(con, format, text, ...)
          {
            con <- file()
            on.exit(close(con))
            writeLines(text, con)
            obj <- import(FileForFormat(con, format), ...)
            obj
          })

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

## Uses XML::parseURI, except first checks for Windows drive letter.
## There are no known URI schemes that are only a single character.
.parseURI <- function(uri) {
  if (.Platform$OS.type == "windows" && grepl("^[A-Za-z]:[/\\]", uri)) {
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

connectionForResource <- function(x, open = "") {
  resource <- decompress(x)
  if (is.character(resource)) {
    if (!nzchar(resource))
      stop("path cannot be an empty string")
    uri <- .parseURI(resource)
    if (uri$scheme != "")
      con <- url(resource)
    else con <- file(resource)
    con <- manage(con)
  } else con <- resource
  if (!isOpen(con) && nzchar(open)) {
    open(con, open)
  }
  con
}

## Connection management (similar to memory management)

manage <- function(con) {
  if (!is.null(attr(con, "finalizerEnv")))
    return(con)
  env <- new.env()
  finalizer <- function(obj) {
    if (exists("con")) {
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
