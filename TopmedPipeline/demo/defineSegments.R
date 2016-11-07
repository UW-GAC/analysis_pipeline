library(dplyr)
library(GenomicRanges)

db <- src_mysql(dbname="hg19", host="genome-mysql.cse.ucsc.edu", user="genome")
chromInfo <- tbl(db, "chromInfo")

dat <- select(chromInfo, chrom, size) %>%
    collect() %>%
    filter(chrom %in% paste0("chr", c(1:22, "X"))) %>%
    mutate(chrom=factor(sub("chr", "", chrom), levels=c(1:22, "X"))) %>%
    arrange(chrom) %>%
    mutate(chrom=as.character(chrom))

gr <- GRanges(seqnames=dat$chrom,
              ranges=IRanges(start=1, end=dat$size))

seg.length <- 1e7
segments <- do.call(c, lapply(gr, function(x) {
    window.start <- seq(start(x), end(x) - seg.length, seg.length)
    window.end <- seq(start(x) + seg.length - 1, end(x), by=seg.length)
    GRanges(seqnames=seqnames(x), IRanges(window.start, window.end))
}))

save(segments, file="data/segments.RData")
