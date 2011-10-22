## code to get twobit support working

library(rtracklayer)
library(Biostrings)

system.time(chr11 <-
            read.DNAStringSet("~/genomes/HuRef/seqs/hs_alt_HuRef_chr11.fa"))

export(chr11, "chr11.2bit")

seqinfo(TwoBitFile("chr11.2bit"))

system.time(chr11_2bit <- import("chr11.2bit"))

identical(chr11_2bit, chr11)

export(subseq(chr11, 1, 100), "test.2bit")
