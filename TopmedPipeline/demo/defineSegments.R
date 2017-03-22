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
chromosomes_hg19 <- gr
save(chromosomes_hg19, file="data/chromosomes_hg19.RData")

library(TopmedPipeline)
seg.length <- 1e7
segments <- defineSegments(seg.length, build="hg19")

# this one will be used for testing
save(segments, file="data/segments.RData")

# this one is the master copy for the pipeline (read by R and python)
writeSegmentFile(segments, "../segments.txt")
