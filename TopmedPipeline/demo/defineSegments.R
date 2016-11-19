library(dplyr)

db <- src_mysql(dbname="hg19", host="genome-mysql.cse.ucsc.edu", user="genome")
chromInfo <- tbl(db, "chromInfo")

dat <- select(chromInfo, chrom, size) %>%
    collect() %>%
    filter(chrom %in% paste0("chr", c(1:22, "X"))) %>%
    mutate(chrom=factor(sub("chr", "", chrom), levels=c(1:22, "X"))) %>%
    arrange(chrom) %>%
    mutate(chrom=as.character(chrom))

library(GenomicRanges)
gr <- GRanges(seqnames=dat$chrom,
              ranges=IRanges(start=1, end=dat$size))

seg.length <- 1e7
segments <- do.call(c, lapply(gr, function(x) {
    window.start <- seq(start(x), end(x), seg.length)
    window.end <- seq(start(x) + seg.length - 1, end(x) + seg.length, by=seg.length)
    GRanges(seqnames=seqnames(x), IRanges(window.start, window.end))
}))

# this one will be used for testing
save(segments, file="data/segments.RData")

seg.df <- as.data.frame(segments) %>%
    dplyr::rename(chromosome=seqnames) %>%
    select(chromosome, start, end)

# this one is the master copy for the pipeline (read by R and python)
write.table(seg.df, file="../segments.txt", quote=FALSE, sep="\t", row.names=FALSE)
