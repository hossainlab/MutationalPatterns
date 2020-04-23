library(GenomicRanges)
library(GenomeInfoDb)


# Drop and Keep Seqlevels 
gr <- GRanges(seqnames = c("chr1", "chr2"), ranges = IRanges(start = 1:2, end = 4:5)) 
seqlevels(gr) <- "chr1"
gr

