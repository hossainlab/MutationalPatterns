library(GenomicRanges)

gr <- GRanges(seqnames = c("chr1", "chr2"), ranges = IRanges(start = 1:2, end = 4:5))
gr
seqinfo(gr)
seqnames(gr)
range(gr)
