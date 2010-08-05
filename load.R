library(RCurl)
library(XML)
library(GenomicRanges)
library(rJava)

files <- c("web.R", "range.R", "trackSet.R", "browser.R", "gff.R", "ucsc.R",
           "bed.R", "wig.R", "io.R")
sapply(files, source)

#track <- import(system.file("inst", "test", "v1.gff", package = "rtracklayer"))
track <- import("../inst/tests/bed.wig")
track@genome <- "hg18"

session <- browserSession("ucsc")
layTrack(session, track)
