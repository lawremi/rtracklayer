library(GGtools)

if (!exists("hmceuB36.2021")) data(hmceuB36.2021)
## condense to founders only
hmFou = hmceuB36.2021[, which(hmceuB36.2021$isFounder)]
## show basic formula fit
f1 = gwSnpTests(genesym("CPNE1")~male, hmFou, chrnum(20))

cpneTrack <- as(f1, "RangedData")
save(cpneTrack, file="../../data/cpneTrack.rda")
