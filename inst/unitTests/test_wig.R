test_wig <- function() {
  test_path <- system.file("tests", package = "rtracklayer")
  test_wig <- file.path(test_path, "step.wig")

  createCorrectRd <- function(si) {
    score <- c(seq(10, 20, by = 2.5), seq(17.5, 10, by = -2.5),
               seq(1000, 100, by = -100))
    start <- c(59104701 + cumsum(c(0, 200, 500, 200, 300, 180, 220, 390, 1180)),
               59108021 + cumsum(c(0, rep(300, 9))))
    width <- c(rep(1, 9), rep(200, 10))
    space <- factor(c(rep("chr19", 9), rep("chr18", 10)), seqlevels(si))
    correct_rd <- RangedData(IRanges(start, width = width), score,
                             space = space)
    if (!any(is.na(genome(si))))
      universe(correct_rd) <- unname(genome(si)[1])
    metadata(ranges(correct_rd))$seqinfo <- si
    correct_rd
  }
  createCorrectUCSC <- function(rd) {  
    track_line <- new("GraphTrackLine", type = "wig", name = "test",
                      description = "test track", visibility = "full",
                      autoScale = FALSE, viewLimits = c(0, 1000),
                      color = c(0L, 200L, 100L),
                      maxHeightPixels = c(100L, 50L, 20L),
                      graphType = "points", priority = 30)
    new("UCSCData", rd, trackLine = track_line)
  }

  correct_rd <- createCorrectRd(Seqinfo(c("chr19", "chr18")))
  correct_ucsc <- createCorrectUCSC(correct_rd)
  
  ## TEST: basic import
  for (asRangedData in c(TRUE, FALSE)) {
    target <- correct_ucsc
    if (!asRangedData)
      target <- as(target, "GRanges")
    test <- import(test_wig, asRangedData = asRangedData)
    checkIdentical(target, test)
    test <- import.wig(test_wig, asRangedData = asRangedData)
    checkIdentical(target, test)
    test_wig_file <- WIGFile(test_wig)
    checkException(import(test_wig_file, format = "bed",
                          asRangedData = asRangedData))
    test <- import(test_wig_file, format = "wig", asRangedData = asRangedData)
    checkIdentical(target, test)
    test <- import(test_wig_file, asRangedData = asRangedData)
    checkIdentical(target, test)
    test_wig_con <- file(test_wig)
    test <- import(test_wig_con, format = "wig", asRangedData = asRangedData)
    checkIdentical(target, test)
    close(test_wig_con)
    test_wig_con <- file(test_wig)
    test <- import(WIGFile(test_wig_con), asRangedData = asRangedData)
    checkIdentical(target, test)
    close(test_wig_con)
  }
  
  ## TEST: 'genome'
  hg19_seqinfo <- SeqinfoForBSGenome("hg19")
  correct_genome <- createCorrectUCSC(createCorrectRd(hg19_seqinfo))
  test <- import(test_wig, genome = "hg19", asRangedData = TRUE)
  checkIdentical(correct_genome, test)
  test <- import(test_wig, genome = "hg19", asRangedData = FALSE)
  checkIdentical(as(correct_genome, "GRanges"), sort(test))

  ## TEST: trackLine = FALSE
  test <- import(test_wig, trackLine = FALSE, asRangedData = TRUE)
  checkIdentical(correct_rd, test)
  test <- import(test_wig, trackLine = FALSE, asRangedData = FALSE)
  checkIdentical(as(correct_rd, "GRanges"), test)

  ## TEST: which
  which <- ranges(correct_rd[3:4,])
  correct_which <- subsetByOverlaps(correct_ucsc, which)
  test <- import(test_wig, which = which, asRangedData = TRUE)
  checkIdentical(correct_which, test)
  test <- import(test_wig, which = which, asRangedData = FALSE)
  checkIdentical(as(correct_which, "GRanges"), test)

  ## TEST: basic export
  test_wig_out <- file.path(tempdir(), "test.wig")
  on.exit(unlink(test_wig_out))
  export(correct_ucsc, test_wig_out)
  test <- import(test_wig_out, asRangedData = TRUE)
  checkIdentical(correct_ucsc, test)
  export.wig(correct_ucsc, test_wig_out)
  test <- import(test_wig_out, asRangedData = TRUE)
  checkIdentical(correct_ucsc, test)
  test_foo_out <- file.path(tempdir(), "test.foo")
  export(correct_ucsc, test_foo_out, format = "wig")
  on.exit(unlink(test_foo_out))
  test <- import(test_wig_out, asRangedData = TRUE)
  checkIdentical(correct_ucsc, test)
  test_wig_out_file <- WIGFile(test_wig_out)
  export(correct_ucsc, test_wig_out_file)
  test <- import(test_wig_out, asRangedData = TRUE)
  checkIdentical(correct_ucsc, test)
  checkException(export(correct_ucsc, test_wig_out_file, format = "gff"))

  ## TEST: append
  correct_ucsc2 <- initialize(correct_ucsc,
                              trackLine = initialize(correct_ucsc@trackLine,
                                name = "test2"))
  export(correct_ucsc2, test_wig_out_file, append = TRUE)
  test <- import(test_wig_out_file, asRangedData = TRUE)
  correct_list <- RangedDataList(test = correct_ucsc,
                                 test2 = correct_ucsc2)
  checkIdentical(correct_list, test)

  ## TEST: track line parameters
  export(correct_ucsc, test_wig_out, name = "test2")
  test <- import(test_wig_out, asRangedData = TRUE)
  checkIdentical(correct_ucsc2, test)

  ## TEST: export trackLine
  export(correct_ucsc, test_wig_out, trackLine = FALSE)
  test <- import(test_wig_out, asRangedData = TRUE)
  checkIdentical(test, correct_rd)

  ## TEST: Plain RangedData / bedGraph
  export.ucsc(correct_rd, test_wig_out) # mixture of steps leads to bedGraph
  test <- import(test_wig_out, asRangedData = TRUE)
  default_line <- new("GraphTrackLine", name = "R Track", type = "bedGraph")
  correct_default <- new("UCSCData", correct_rd, trackLine = default_line)
  checkIdentical(test, correct_default)
  export.ucsc(correct_rd[2], test_wig_out)
  test <- import(test_wig_out, asRangedData = TRUE)
  default_line <- new("GraphTrackLine", name = "R Track", type = "wig")
  correct_default <- new("UCSCData", correct_rd[2], trackLine = default_line)
  seqinfo(correct_default) <- seqinfo(correct_default)["chr18"]
  checkIdentical(test, correct_default)
  
  ## TEST: RangedDataList
  export(correct_list, test_wig_out)
  test <- import(test_wig_out, asRangedData = TRUE)
  checkIdentical(correct_list, test)
  
  ## TEST: gzip
  test_wig_gz <- paste(test_wig_out, ".gz", sep = "")
  on.exit(unlink(test_wig_gz))
  export(correct_ucsc, test_wig_gz)
  test <- import(test_wig_gz, asRangedData = TRUE)
  checkIdentical(correct_ucsc, test)
  export(correct_ucsc2, test_wig_gz, append = TRUE)
  test <- import(test_wig_gz, asRangedData = TRUE)
  checkIdentical(correct_list, test)
  
  ## TEST: Using connection to add comment header
  test_wig_con <- file(test_wig_out)
  open(test_wig_con, "w")
  comment <- "# test comment"
  writeLines(comment, test_wig_con)
  export(correct_ucsc, test_wig_con)
  close(test_wig_con)
  checkIdentical(comment, readLines(test_wig_out, n = 1))
  test <- import(test_wig_out, asRangedData = TRUE)
  checkIdentical(correct_ucsc, test)
}
