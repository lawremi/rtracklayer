library(rtracklayer)
files <- dir(system.file("tests", package = "rtracklayer"), pattern = "[^R~]$",
             full.names=TRUE)
testExport <- function(file) {
  track <- import(file)
  export(track, format = "ucsc")
}
options(error=recover)
sapply(files, testExport)

## test GRanges imports
sapply(files[-8], import, asRangedData = FALSE)
