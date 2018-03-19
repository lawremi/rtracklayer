###################################################
### chunk number 1: rtl-init
###################################################
library(rtracklayer)
data(targets)

###################################################
### chunk number 2: rtl-miRNA-track
###################################################
targetTrack <- makeGRangesFromDataFrame(targets,
                   keep.extra.columns=TRUE)


###################################################
### chunk number 3: rtl-export eval=FALSE
###################################################
## export(targetTrack, "targets.wig")


###################################################
### chunk number 4: rtl-ucsc-start
###################################################
session <- browserSession()
genome(session) <- "hg18"


###################################################
### chunk number 5: rtl-ucsc-lay
###################################################
session$targets <- targetTrack


###################################################
### chunk number 6: rtl-ucsc-view eval=FALSE
###################################################
top <- targetTrack$target == targets$target[1]
range <- targetTrack[top,] * -10
view <- browserView(session, range,
                    hide = c("refGene", "mgcFullMrna", "intronEst"),
                    dense = "knownGene", squish = "cons44way")

