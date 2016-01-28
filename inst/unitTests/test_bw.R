test_bw <- function() {
  if (.Platform$OS.type == "windows")
    return()
  
  test_path <- system.file("tests", package = "rtracklayer")

  test_bw <- file.path(test_path, "test.bw")

  ir <- as(PartitioningByWidth(rep(300, 9)), "IRanges")
  space <- factor(c(rep("chr2", 5), rep("chr19", 4)), c("chr2", "chr19"))
  score <- seq(-1, 1, length = 9)
  correct_fixed <- GRanges(space, ir, score = score)
  si <- SeqinfoForBSGenome("hg19")
  seqlengths(correct_fixed) <- seqlengths(si)[levels(space)]

  test <- import(test_bw)
  checkIdentical(test, correct_fixed)

  test_bw_out <- file.path(tempdir(), "test_out.bw")
  export(correct_fixed, test_bw_out)
  on.exit(unlink(test_bw_out))
  test <- import(test_bw_out)
  checkIdentical(test, correct_fixed)

  export.bw(correct_fixed, test_bw_out)
  test <- import.bw(test_bw_out)
  checkIdentical(test, correct_fixed)
  
  correct_bedgraph <- correct_fixed
  width(correct_bedgraph) <- seq(1, 300, length = 9)
  export(correct_bedgraph, test_bw_out)
  test <- import(test_bw_out)
  checkIdentical(test, correct_bedgraph)
  
  ## TEST: 'which'
  which <- GRanges(c("chr2", "chr2"), IRanges(c(1, 300), c(400, 1000)))
  correct_which <- subsetByOverlaps(correct_bedgraph, which)
  ranges(correct_which) <- ranges(intersect(correct_which, which))
  test <- import(test_bw_out, which = which)
  checkIdentical(test, correct_which)

  ## TEST: BigWigSelection (range, no score)
  test <- import(test_bw_out,
                 selection = BigWigSelection(which, colnames = character()))
  correct_which <- correct_which[, character()]
  checkIdentical(test, correct_which)
  
  ## TEST: empty which
  which <- GRanges()
  correct_which <- subsetByOverlaps(correct_bedgraph, which)
  test <- import(test_bw_out, which = which)
  checkIdentical(test, correct_which)  
  
  ## TEST: non-UCSC naming
  correct_ncbi <- correct_bedgraph
  seqlevels(correct_ncbi) <- sub("chr", "", seqlevels(correct_ncbi))
  export(correct_ncbi, test_bw_out)
  test <- import(test_bw_out)
  checkIdentical(test, correct_ncbi)

  ## TEST: as="RleList"
  correct_cov <- coverage(correct_ncbi, weight="score")
  test <- import(test_bw_out, as="RleList")
  checkIdentical(correct_cov, test)

  ## TEST: export RleList
  export(correct_cov, test_bw_out)
  test <- import(test_bw_out, as="RleList")
  checkIdentical(correct_cov, test)

  ## TEST: export/import NumericList
  correct_cov_short <- correct_cov[correct_cov != 0L]
  correct_int <- as(correct_cov_short, "NumericList")
  which <- GRanges(names(correct_int), IRanges(1, elementNROWS(correct_int)))
  metadata(correct_int) <- list(ranges=which)
  export(correct_int, test_bw_out)

  test <- import(test_bw_out, as="NumericList")
  checkIdentical(correct_int, test)
  test <- import(test_bw_out, which=which[1], as="NumericList")
  checkIdentical(elementNROWS(correct_int[1]), elementNROWS(test))
  test <- import(test_bw_out, which=which[1:2], as="NumericList")
  checkIdentical(correct_int, test)
}
