# import/export

setGeneric("export",
           function(object, con, format, ...) standardGeneric("export"))

setMethod("export", c(format = "character"),
          function(object, con, format, ...)
          {
            fun <- try(match.fun(paste("export", format, sep=".")), TRUE)
            if (is.character(fun))
              stop("No export function for '", format, "' found")
            wasOpen <- TRUE
            if (!isOpen(con)) {
              open(con, "w")
              wasOpen <- FALSE
            }
            fun(object, con, ...)
            if (!wasOpen)
              close(con)
          })
setMethod("export", c(con = "missing", format = "character"),
          function(object, con, format, ...)
          {
            con <- file()
            export(object, con, format, ...)
            text <- readLines(con, warn = FALSE)
            close(con)
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
            con <- file(con)
            export(object, con, format, ...)
          })

setGeneric("import",
           function(con, format, text, ...) standardGeneric("import"))

setMethod("import", c(format = "character"),
          function(con, format, text, ...)
          {
            fun <- try(match.fun(paste("import", format, sep=".")), TRUE)
            if (is.character(fun))
              stop("No import function for '", format, "' found")
            wasOpen <- TRUE
            if (!isOpen(con)) {
              open(con, "r")
              wasOpen <- FALSE
            }
            object <- fun(con, ...)
            if (!wasOpen)
              close(con)
            object
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
            con <- file(con)
            import(con, format, ...)
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

# utilities
file_ext <- function(con) gsub(".*\\.([^.]*)$", "\\1", con)
