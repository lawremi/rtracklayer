test_bedGraph <- function() {
  test_path <- system.file("tests", package = "rtracklayer")
  test_bg <- file.path(test_path, "test.bedGraph")

  createCorrectGR <- function(si) {
    part <- PartitioningByWidth(rep(300, 9))
    ir <- shift(IRanges(start(part), end(part)), 59102000)
    score <- seq(-1, 1, by = 0.25)
    space <- factor(c(rep("chr19", 6), "chr17", rep("chr18", 2)), seqlevels(si))
    correct_gr <- GRanges(space, ir, score = score)
    if (!any(is.na(genome(si))))
      genome(correct_gr) <- unname(genome(si)[1])
    seqinfo(correct_gr) <- si
    correct_gr
  }
  createCorrectUCSC <- function(gr) {  
    track_line <- new("GraphTrackLine", type = "bedGraph",
                      name = "bedGraph track",
                      description = "Test", visibility = "full",
                      alwaysZero = TRUE, gridDefault = FALSE, yLineMark = 10,
                      windowingFunction = "mean",
                      color = c(200L, 100L, 0L),
                      altColor = c(0L, 100L, 200L), priority = 20)
    new("UCSCData", gr, trackLine = track_line)
  }

  correct_gr <- createCorrectGR(Seqinfo(c("chr19", "chr17", "chr18")))
  correct_ucsc <- createCorrectUCSC(correct_gr)
  
  ## TEST: basic import
  test <- import(test_bg)
  checkIdentical(correct_ucsc, test)
  test <- import.bedGraph(test_bg)
  checkIdentical(correct_ucsc, test)
  test_bg_file <- BEDGraphFile(test_bg)
  checkException(import(test_bg_file, format = "wig"))
  test <- import(test_bg_file, format = "bedGraph")
  checkIdentical(correct_ucsc, test)
  test <- import(test_bg_file)
  checkIdentical(correct_ucsc, test)
  test_bg_con <- file(test_bg)
  test <- import(test_bg_con, format = "wig")
  checkIdentical(correct_ucsc, test)
  test_bg_con <- file(test_bg)
  test <- import(WIGFile(test_bg_con))
  checkIdentical(correct_ucsc, test)

  ## TEST: 'genome'
  hg19_seqinfo <- SeqinfoForBSGenome("hg19")
  correct_genome <- createCorrectUCSC(createCorrectGR(hg19_seqinfo))
  test <- import(test_bg, genome = "hg19")
  checkIdentical(correct_genome, test)

  ## TEST: trackLine = FALSE
  test <- import(test_bg, trackLine = FALSE)
  checkIdentical(correct_gr, test)

  ## TEST: which
  correct_which <- subsetByOverlaps(correct_ucsc, correct_gr[3:4])
  test <- import(test_bg, which = correct_which)
  checkIdentical(correct_which, test)

  ## TEST: basic export
  test_bg_out <- file.path(tempdir(), "test.bedGraph")
  on.exit(unlink(test_bg_out))
  export(correct_ucsc, test_bg_out)
  test <- import(test_bg_out)
  checkIdentical(correct_ucsc, test)
  export.bedGraph(correct_ucsc, test_bg_out)
  test <- import(test_bg_out)
  checkIdentical(correct_ucsc, test)
  test_foo_out <- file.path(tempdir(), "test.foo")
  export(correct_ucsc, test_foo_out, format = "bedGraph")
  on.exit(unlink(test_foo_out))
  test <- import(test_bg_out)
  checkIdentical(correct_ucsc, test)
  test_bg_out_file <- BEDGraphFile(test_bg_out)
  export(correct_ucsc, test_bg_out_file)
  test <- import(test_bg_out)
  checkIdentical(correct_ucsc, test)
  checkException(export(correct_ucsc, test_bg_out_file, format = "gff"))

  ## TEST: append
  correct_ucsc2 <- initialize(correct_ucsc,
                              trackLine = initialize(correct_ucsc@trackLine,
                                                     name = "test2"))
  export(correct_ucsc2, test_bg_out_file, append = TRUE)
  test <- import(test_bg_out_file)
  correct_list <- GRangesList("bedGraph track" = correct_ucsc,
                              test2 = correct_ucsc2,
                              compress = FALSE)
  checkIdentical(correct_list, test)

  ## TEST: track line parameters
  export(correct_ucsc, test_bg_out, name = "test2")
  test <- import(test_bg_out)
  checkIdentical(correct_ucsc2, test)

  ## TEST: export trackLine
  export(correct_ucsc, test_bg_out, trackLine = FALSE)
  test <- import(test_bg_out)
  checkIdentical(test, correct_gr)
  
  ## TEST: SimpleGRangesList
  export(correct_list, test_bg_out)
  test <- import(test_bg_out)
  checkIdentical(correct_list, test)
  
  ## TEST: gzip
  test_bg_gz <- paste(test_bg_out, ".gz", sep = "")
  on.exit(unlink(test_bg_gz))
  export(correct_ucsc, test_bg_gz)
  test <- import(test_bg_gz)
  checkIdentical(correct_ucsc, test)
  export(correct_ucsc2, test_bg_gz, append = TRUE)
  test <- import(test_bg_gz)
  checkIdentical(correct_list, test)
  
  ## TEST: Using connection to add comment header
  test_bg_con <- file(test_bg_out)
  open(test_bg_con, "w")
  comment <- "# test comment"
  writeLines(comment, test_bg_con)
  export(correct_ucsc, test_bg_con)
  close(test_bg_con)
  checkIdentical(comment, readLines(test_bg_out, n = 1))
  test <- import(test_bg_out)
  checkIdentical(correct_ucsc, test)
}
