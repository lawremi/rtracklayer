### =========================================================================
### BAM support (wrappers around Rsamtools)
### -------------------------------------------------------------------------

setMethod("import", "BamFile",
          function(con, format, text, use.names = FALSE,
                   param = ScanBamParam(...), ...)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            readGAlignmentsFromBam(con, use.names = use.names, param = param)
          })

setMethod("export", c("GAlignments", "BamFile"),
          function(object, con, format, index = TRUE) {
            sam_path <- paste(file_path_sans_ext(path(con)), ".sam", sep = "")
            sam_con <- file(sam_path, "w")
            on.exit(close(sam_con))
            si <- seqinfo(object)
            has_info <-
              seqlevels(si)[!is.na(seqlevels(si)) & !is.na(seqlengths(si))]
            si <- si[has_info]
            if (length(si)) {
              header <- paste0("@SQ",
                               "\tSN:", seqlevels(si),
                               "\tLN:", seqlengths(si))
              has_genome <- !is.na(genome(si))
              header[has_genome] <-  paste0(header[has_genome], "\tAS:",
                                            genome(si)[has_genome])
              writeLines(header, sam_con)
            }
            emd <- mcols(object)
            aln <- paste(if (!is.null(names(object))) names(object) else "*",
                         if (!is.null(emd[["flag"]])) emd[["flag"]] else
                           ifelse(strand(object) == "-", "16", "0"),
                         seqnames(object), start(object),
                         if (!is.null(emd[["mapq"]])) emd[["mapq"]] else "255",
                         cigar(object),
                         if (!is.null(emd[["mrnm"]])) emd[["mrnm"]] else "*",
                         if (!is.null(emd[["mpos"]])) emd[["mpos"]] else "0",
                         if (!is.null(emd[["isize"]])) emd[["isize"]] else "0",
                         if (is(object, "GappedReads")) object@qseq
                           else if (!is.null(emd[["seq"]])) emd[["seq"]]
                           else "*",
                         if (!is.null(emd[["qual"]])) emd[["qual"]] else "*",
                         sep = "\t")
            custom <- emd[nchar(names(emd)) == 2L]
            if (length(custom) > 0L) {
              type.map <- c(integer = "i", numeric = "f", character = "Z",
                            factor = "Z")
              custom.class <- sapply(custom, class)
              custom.type <- type.map[custom.class]
              unknown.class <- custom.class[is.na(custom.type)]
              if (length(unknown.class) > 0L) {
                warning("these classes are not yet valid for BAM tag export: ",
                        paste(unknown.class, collapse=", "))
                custom <- custom[!is.na(custom.type)]
              }
              tags <- mapply(paste0, names(custom), ":", custom.type, ":",
                             as.list(custom), SIMPLIFY=FALSE)
              aln <- do.call(paste, c(list(aln), tags, sep = "\t"))
            }
            writeLines(aln, sam_con)
            close(sam_con)
            on.exit()
            bam <- asBam(sam_path, file_path_sans_ext(sam_path),
                         overwrite = TRUE, indexDestination = index)
            unlink(sam_path)
            invisible(bam)
          })

setMethod("export", c("ANY", "BamFile"),
          function(object, con, format, ...) {
            export(as(object, "GAlignments"), con, ...)
          })

