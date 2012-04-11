test_bw <- function() {
  if (.Platform$OS.type == "windows")
    return()
  
  test_path <- system.file("tests", package = "rtracklayer")

  test_bw <- file.path(test_path, "test.bw")

  ir <- as(PartitioningByWidth(rep(300, 9)), "IRanges")
  space <- factor(c(rep("chr2", 5), rep("chr19", 4)), c("chr2", "chr19"))
  score <- seq(-1, 1, length = 9)
  correct_fixed <- RangedData(ir, score, space = space)
  ## convert the Compressed lists to Simple equivalents for identity with C code
  ranges(correct_fixed) <- RangesList(as.list(ranges(correct_fixed)))
  values(correct_fixed) <- SplitDataFrameList(as.list(values(correct_fixed)),
                                              compress = FALSE)
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
  ranges(correct_which) <- intersect(ranges(correct_which),
                                     as(which, "RangesList"))
  test <- import(test_bw_out, which = which)
  checkIdentical(test, correct_which)

  ## TEST: BigWigSelection (range, no score)
  test <- import(test_bw_out,
                 selection = BigWigSelection(which, colnames = character()))
  correct_which <- correct_which[, character()]
  checkIdentical(test, correct_which)
  
  ## TEST: empty which
  which <- RangesList()
  correct_which <- subsetByOverlaps(correct_bedgraph, which)
  test <- import(test_bw_out, which = which)
  checkIdentical(test, correct_which)  
  
  ## TEST: asRangedData=FALSE
  test <- import(test_bw, asRangedData = FALSE)
  checkIdentical(test, as(correct_fixed, "GRanges"))
  
  ## TEST: non-UCSC naming
  correct_ncbi <- correct_bedgraph
  seqlevels(correct_ncbi) <- sub("chr", "", seqlevels(correct_ncbi))
  names(correct_ncbi) <- sub("chr", "", names(correct_ncbi))
  export(correct_ncbi, test_bw_out)
  test <- import(test_bw_out)
  checkIdentical(test, correct_ncbi)  
}
