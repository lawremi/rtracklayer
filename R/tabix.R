### =========================================================================
### Tabix support
### -------------------------------------------------------------------------
###

setMethod("import", c("TabixFile", "character"),
          function(con, format, text,
                   which = as(seqinfoForGenome(genome), "GenomicRanges"),
                   genome = NA, ...)
          {
            if (missing(which) && is.na(genome))
              stop("'which' or 'genome' must be specified")
            lines <- scanTabix(con, ..., which)
            buffer <- file()
            writeLines(lines, buffer)
            fun <- .importForFormat(format)
            if (is.null(fun)) {
              header <- headerTabix(con)
              do.call(import.table,
                      c(buffer, genome = genome, header$indexColumns,
                        skip = header$skip, list(...)))
            }
            fun(buffer, genome = genome, ...)
          })

setMethod("import.gff", "TabixFile",
          function(con, ...)
          {
            import(con, format = "gff", ...)
          })

setMethod("import.bed", "TabixFile",
          function(con, ...)
          {
            import(con, format = "bed", ...)
          })

setMethod("import.bedGraph", "TabixFile",
          function(con, ...)
          {
            import(con, format = "bedGraph", ...)
          })

setMethod("import.bed15", "TabixFile",
          function(con, ...)
          {
            import(con, format = "bed15", ...)
          })
