library(RCurl)
library(XML)
library(Biobase)
library(rJava)

files <- c("web.R", "segment.R", "trackSet.R", "browser.R", "gff.R", "ucsc.R",
           "bed.R", "wig.R", "io.R", "argo.R")
sapply(files, source)

#track <- import(system.file("inst", "test", "v1.gff", package = "rtracklayer"))
track <- import("../inst/tests/bed.wig")
track@genome <- "hg18"

#session <- browseGenome(browser = "argo")
session <- browserSession("argo")
layTrack(session, track)
