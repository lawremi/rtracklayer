library(rtracklayer)
library(RUnit)

test_wig <- function() {
  test_path <- system.file("tests", package = "rtracklayer")
  test_wig <- file.path(test_path, "step.wig")

  createCorrectRd <- function(si) {
    score <- c(seq(10, 20, by = 2.5), seq(17.5, 10, by = -2.5),
               seq(1000, 100, by = -100))
    start <- c(59304701 + cumsum(c(0, 200, 500, 200, 300, 180, 220, 690, 1180)),
               59308021 + cumsum(c(0, rep(300, 9))))
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
  test <- import(test_wig)
  checkIdentical(test, correct_ucsc)
  test <- import.wig(test_wig)
  checkIdentical(test, correct_ucsc)
  test_wig_file <- WIGFile(test_wig)
  test <- import(test_wig_file)
  checkIdentical(import(test_wig_file, format = "wig"), correct_ucsc)
  checkException(import(test_wig_file, format = "bed"))
  checkIdentical(test, correct_ucsc)
  test_wig_con <- file(test_wig)
  test <- import(test_wig_con, format = "wig")
  checkIdentical(test, correct_ucsc)
  close(test_wig_con)
  test_wig_con <- file(test_wig)
  test <- import(BEDFile(test_wig_con))
  checkIdentical(test, correct_ucsc)
  close(test_wig_con)
  
  ## TEST: 'genome'
  hg19_seqinfo <- SeqinfoForBSGenome("hg19")
  correct_genome <- createCorrectUCSC(createCorrectRd(hg19_seqinfo))
  test <- import(test_wig, genome = "hg19")
  checkIdentical(correct_genome, test)

  ## TEST: trackLine = FALSE
  ## test <- import(test_wig, trackLine = FALSE)
  ## checkIdentical(test, correct_rd)

  ## TEST: asRangedData = FALSE
  ## correct_gr <- as(correct_ucsc, "GRanges")
  ## test <- import(test_wig, asRangedData = FALSE)
  ## checkIdentical(correct_gr, sort(test))
}
