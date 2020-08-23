test_bb <- function() {
  if (.Platform$OS.type == "windows")
    return()

  test_path <- system.file("tests", package = "rtracklayer")
  test_bb <- file.path(test_path, "test.bb")
  start <- c(237640, 521500 ,565725, 565900, 566760,
             119905, 122525, 173925, 179865, 180185)
  ir <- IRanges(start, width = 151)
  space <- factor(c(rep("chr1", 5), rep("chr10", 5)))
  name <- rep(".", 10)
  score <- seq.int(70L, 700L, length = 10)
  signalValue <- seq(10, 100, length = 10)
  peak <- rep(-1L, 10)
  correct_fixed <- GRanges(space, ir, name = name, score = score,
                           signalValue = signalValue , peak = peak)
  si <- SeqinfoForBSGenome("hg19")
  seqlengths(correct_fixed) <- seqlengths(si)[levels(space)]

  ## TEST: import whole file
  test <- import(test_bb)
  checkIdentical(test, correct_fixed)

  ## TEST: 'which'
  which <- GRanges(c("chr10"), IRanges(c(180185, 180335)))
  correct_which <- subsetByOverlaps(correct_fixed, which)
  test <- import(test_bb, which = which)
  checkIdentical(test, correct_which)

  ## TEST: empty which
  which <- GRanges()
  correct_which <- subsetByOverlaps(correct_fixed, which)
  test <- import(test_bb, which = which)
  checkIdentical(test, correct_which)

  ## TEST: BigBedSelection (GRanges, no field)
  which <- GRanges(c("chr10"), IRanges(c(180185, 180335)))
  test <- import(test_bb,
                 selection = BigBedSelection(which, colnames = character()))
  correct_subset <- subsetByOverlaps(correct_fixed, which)
  correct_which <- correct_subset[, character()]
  correct_which@elementMetadata <- DataFrame()
  checkIdentical(test, correct_which)

  ## TEST: BigBedSelection (GRanges, 1 default field)
  test <- import(test_bb,
                 selection = BigBedSelection(which, colnames = c("name")))
  correct_which <- correct_subset[, c("name")]
  checkIdentical(test, correct_which)

  ## TEST: BigBedSelection (GRanges, 1 extra field)
  test <- import(test_bb,
                 selection = BigBedSelection(which, colnames = c("peak")))
  correct_which <- correct_subset[, c("peak")]
  checkIdentical(test, correct_which)

  ## TEST: BigBedSelection (GRanges, 1 default field and 1 extra field)
  colnames <- c("name", "peak")
  test <- import(test_bb,
                 selection = BigBedSelection(which, colnames =colnames))
  correct_which <- correct_subset[, colnames]
  checkIdentical(test, correct_which)

  # TEST: export
  test_bb_out <- file.path(tempdir(), "test_out.bb")
  export(correct_fixed, test_bb_out)
  on.exit(unlink(test_bb_out))
  test <- import(test_bb_out)
  checkIdentical(test, correct_fixed)
}
