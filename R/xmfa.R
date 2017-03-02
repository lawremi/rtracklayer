### =========================================================================
### Input/output of XMFA files
### -------------------------------------------------------------------------

setClass("XMFAFile", contains = "RTLFile")
XMFAFile <- function(resource) {
    new("XMFAFile", resource = resource)
}

setMethod("import", "XMFAFile", function(con, format, text) {
    lines <- readLines(resource(con))
    lines <- lines[!grepl("^#", lines)]
        header.pos <- grep("^>", lines)
    header.df <-
        strcapture("^> .*?:([[:digit:]]+)-([[:digit:]]+) ([+-]) (.*)",
                   lines[header.pos],
                   DataFrame(start=integer(), end=integer(),
                             strand=character(),
                             seqnames=character()))
    
    lcp.pos <- grep("^=", lines)
    seq.breaks <- sort(c(header.pos, lcp.pos))
    seq.ranges <- IRanges(head(seq.breaks+1L, -1L),
                          tail(seq.breaks-1L, -1L))
    seq.ranges <- seq.ranges[width(seq.ranges) > 0L]
    seqs <-
        DNAStringSet(unstrsplit(extractList(lines, seq.ranges), ""))
    dels <- deletionsFromGaps(seqs)
    unaligned.lcp <- seqs[gaps(dels, 1L, width(seqs))]
    
    range <- with(header.df, IRanges(start, end))
    ord.lcp <- order(start(range))
    unaligned.seqname <- split(unaligned.lcp[ord.lcp],
                               header.df$seqnames[ord.lcp])
    unaligned <- unstrsplit(unaligned.seqname, "")[header.df$seqnames]
    
    lcp.breaks <- c(1L, lcp.pos)
    lcp <- findInterval(header.pos, lcp.breaks)

    mismatch <- mismatchesFromLCPDiscordance(seqs, lcp)
    mismatch <- mapMismatches(mismatch, range, dels)
    
    inverted <- header.df$strand == "-"
    
    aln <- AlignedXStringSet(unaligned,
                             range,
                             mismatch,
                             dels,
                             inverted)
    mcols(aln)$seqname <- header.df$seqname
    unname(splitAsList(aln, lcp))
})


mismatchesFromLCPDiscordance <- function(seqs, lcp) {
    seqs.lcp <- split(seqs, lcp)
    mm <- IntegerList(lapply(seqs.lcp, discordantPositions))
    rep(mm, lengths(seqs.lcp))
}

discordantPositions <- function(x) {
    mat <- consensusMatrix(x)
    baseMat <- mat[DNA_BASES,]
    which(colSums(baseMat == length(x) - mat["-",][col(baseMat)]) == 0L)
}
