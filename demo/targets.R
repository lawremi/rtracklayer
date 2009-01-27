###################################################
### chunk number 1: rtl-init
###################################################
library(rtracklayer)
data(targets)

###################################################
### chunk number 2: rtl-miRNA-track
###################################################
targetTrack <- with(targets, 
    RangedData(IRanges(start, end), 
               target, strand, space = chrom))


###################################################
### chunk number 3: rtl-export eval=FALSE
###################################################
## export(targetTrack, "targets.wig")


###################################################
### chunk number 4: rtl-ucsc-start
###################################################
session <- browserSession()


###################################################
### chunk number 5: rtl-ucsc-lay
###################################################
session$targets <- targetTrack


###################################################
### chunk number 6: rtl-ucsc-view eval=FALSE
###################################################
range <- ranges(targetTrack[1,]) * -10
view <- browserView(session, range,
                    pack = "targets")


