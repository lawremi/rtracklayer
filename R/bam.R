### =========================================================================
### BAM support (wrappers around Rsamtools)
### -------------------------------------------------------------------------

setMethod("import", "BamFile",
          function(con, format, text, use.names = FALSE,
                   param = ScanBamParam(...), ...)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            readBamGappedAlignments(con, use.names = use.names, param = param)
          })

setMethod("export", c("GappedAlignments", "BamFile"),
          function(object, con, format) {
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
            emd <- values(object)
            aln <- paste(if (!is.null(names(object))) names(object) else "*",
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
            writeLines(aln, sam_con)
            close(sam_con)
            on.exit()
            bam <- asBam(sam_path, file_path_sans_ext(sam_path))
            unlink(sam_path)
            invisible(bam)
          })

setMethod("export", c("ANY", "BamFile"),
          function(object, con, format, ...) {
            export(as(object, "GappedAlignments"), con, ...)
          })

