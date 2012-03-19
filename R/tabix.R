### =========================================================================
### Tabix support
### -------------------------------------------------------------------------
###

setMethod("import", c("TabixFile", "character"),
          function(con, format, text,
                   which = as(seqinfoForGenome(genome), "GenomicRanges"),
                   genome = NA, header = TRUE, ...)
          {
            if (missing(which) && is.na(genome))
              stop("'which' or 'genome' must be specified")
            buffer <- queryForResource(con, which, header = header)
            file <- try(fileForFormat(buffer, format), silent = TRUE)
            if (is(file, "try-error")) {
              tabixHeader <- headerTabix(con)
              args <- list(...)
              if (!header && is.null(args$skip))
                args$skip <- tabixHeader$skip
              if (is.null(args$comment.char))
                args$comment.char <- tabixHeader$comment
              do.call(import.trackTable,
                      c(buffer, genome = genome, tabixHeader$indexColumns,
                        args))
            } else import(file, genome = genome, ...)
          })

setMethod("import", c("TabixFile", "missing"),
          function(con, format, text,
                   which = as(seqinfoForGenome(genome), "GenomicRanges"),
                   genome = NA, header = TRUE, ...)
          {
            format <- file_ext(file_path_sans_ext(path(con)))
            callGeneric()
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
