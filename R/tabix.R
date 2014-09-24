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
              skip <- tabixHeader$skip
              if (header) {
                skip <- skip - 1L
              }
              do.call(import.tabSeparated,
                      c(list(buffer), genome = genome, tabixHeader$indexColumns,
                        skip = skip, header = header, args))
            } else import(file, genome = genome, ...)
          })

setMethod("import", c("TabixFile", "missing"),
          function(con, format, text, ...)
          {
            import(con, file_ext(file_path_sans_ext(path(con))), ...)
          })

setGeneric("exportToTabix",
           function(object, con, ...) standardGeneric("exportToTabix"))

setMethod("exportToTabix", c("ANY", "character"),
          function(object, con, ...) {
            con <- TabSeparatedFile(con)
            export(sort(object, ignore.strand=TRUE), con,
                   row.names=FALSE, col.names=TRUE, ...)
            indexTrack(con, seq=1L, start=2L, end=3L, skip=1L)
          })
