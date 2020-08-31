### =========================================================================
### BAM support (wrappers around Rsamtools)
### -------------------------------------------------------------------------

setMethod("import", "BamFile",
          function(con, format, text, paired = FALSE, use.names = FALSE,
                   param = ScanBamParam(...), genome = NA_character_, ...)
          {
            if (!missing(format))
                checkArgFormat(con, format)
            stopifnot(isSingleStringOrNA(genome) || is(genome, "Seqinfo"))
            stopifnot(isTRUEorFALSE(paired))
            if (paired) {
                ans <- readGAlignmentPairs(con, use.names = use.names,
                                           param = param)
            } else {
                ans <- readGAlignments(con, use.names = use.names,
                                       param = param)
            }
            if (isSingleStringOrNA(genome)) {
                genome <- Seqinfo(genome=genome)
            }
            seqinfo(ans) <- merge(seqinfo(ans), genome)
            ans
          })

fillColumn <- function(x, filler) {
    if (is.null(x))
        filler
    else if (anyNA(x))
        ifelse(is.na(x), filler, x)
    else x
}

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
            aln <- paste(fillColumn(names(object), "*"),
                         fillColumn(emd[["flag"]],
                                    ifelse(strand(object) == "-", "16", "0")),
                         seqnames(object), start(object),
                         fillColumn(emd[["mapq"]], "255"),
                         cigar(object),
                         fillColumn(emd[["mrnm"]], "*"),
                         fillColumn(emd[["mpos"]], "0"),
                         fillColumn(emd[["isize"]], "0"),
                         if (is(object, "GappedReads")) object@qseq
                         else fillColumn(emd[["seq"]], "*"),
                         fillColumn(emd[["qual"]], "*"),
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

setMethod("export", c("GAlignmentPairs", "BamFile"),
          function(object, con, format, ...) {
            ga <- as(object, "GAlignments")
            collate_subscript <-
              S4Vectors:::make_XYZxyz_to_XxYyZz_subscript(length(object))
            getMateAttribute <- function(FUN) {
              c(FUN(last(object)), FUN(first(object)))[collate_subscript]
            }
            if (is.null(mcols(ga)$mrnm)) {
              mcols(ga)$mrnm <- getMateAttribute(seqnames)
            }
            if (is.null(mcols(ga)$mpos)) {
              mcols(ga)$mpos <- getMateAttribute(start)
            }
### FIXME: we cannot infer whether the pair is 'proper' (0x2)
### nor whether the alignment is 'primary' (0x100)
            if (is.null(mcols(ga)$flag)) {
              mcols(ga)$flag <- 0x1 + # all paired
                ifelse(strand(ga) == "-", 0x10, 0) +
                  ifelse(getMateAttribute(strand) == "-", 0x20, 0) +
                    c(0x40, 0x80) # left vs. right
            }
            if (is.null(names(ga))) {
              names(ga) <- as.character(rep(seq_len(length(object)), each=2L))
            }
            export(ga, con, ...)
          })

setMethod("export", c("ANY", "BamFile"),
          function(object, con, format, ...) {
            export(as(object, "GAlignments"), con, ...)
          })

