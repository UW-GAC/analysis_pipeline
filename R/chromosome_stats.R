library(dplyr)

db <- src_mysql(dbname="hg19", host="genome-mysql.cse.ucsc.edu", user="genome")
chromInfo <- tbl(db, "chromInfo")

dat <- select(chromInfo, chrom, size) %>%
    collect() %>%
    filter(chrom %in% paste0("chr", c(1:22, "X"))) %>%
    mutate(chrom=sub("chr", "", chrom))
