### =========================================================================
### Tabix support
### -------------------------------------------------------------------------
###

setMethod("import", c("TabixFile", "character"),
          function(con, format, text,
                   which = if (is.na(genome)) NULL
                           else as(seqinfoForGenome(genome), "GenomicRanges"),
                   genome = NA, header = TRUE, ...)
          {
            buffer <- queryForResource(con, which, header = header)
            on.exit(release(buffer))
            file <- try(FileForFormat(buffer, format), silent = TRUE)
            if (is(file, "try-error")) {
              tabixHeader <- headerTabix(con)
              args <- list(...)
              if (is.null(args$comment.char))
                args$comment.char <- tabixHeader$comment
              do.call(import.trackTable,
                      c(list(buffer), genome = genome, tabixHeader$indexColumns,
                        args))
            } else import(file, genome = genome, ...)
          })

setMethod("import", c("TabixFile", "missing"),
          function(con, format, text, ...)
          {
            import(con, file_ext(file_path_sans_ext(path(con))), ...)
          })
