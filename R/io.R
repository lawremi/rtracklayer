## import/export

## Need to dispatch on S3 "connection"

##setOldClass("connection")

.connectionClasses <- c("file", "url", "gzfile", "bzfile", "unz", "pipe",
                        "fifo", "sockconn", "terminal", "textConnection")
apply(cbind(.connectionClasses, "connection"), 1, setOldClass,
      where = environment())
setClassUnion("characterORconnection", c("character", "connection"))

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
            fun <- .importForFormat(format)
            fun(con, ...)
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
setGeneric("import.gz",
           function(con, ...) standardGeneric("import.gz"))

setMethod("import.gz", "character", function(con, ...) {
  import(gzfile(con), file_ext(sub("\\.gz$", "", con)), ...)
})

setMethod("import.gz", "connection", function(con, ...) {
  import(gzcon(con), ...)
})

## utilities
file_ext <- function(con) gsub(".*\\.([^.]*)$", "\\1", con)
