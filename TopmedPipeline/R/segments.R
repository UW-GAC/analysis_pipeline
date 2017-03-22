defineSegments <- function(seg.length=1e7, build="hg19") {
    gr <- get(data(list=paste("chromosomes", build, sep="_"),
                   package="TopmedPipeline", envir=environment()))
    do.call(c, lapply(gr, function(x) {
        window.start <- seq(start(x), end(x), seg.length)
        window.end <- seq(start(x) + seg.length - 1, end(x) + seg.length, by=seg.length)
        GRanges(seqnames=seqnames(x), IRanges(window.start, window.end))
    }))
}

writeSegmentFile <- function(segments, file) {
    seg.df <- as.data.frame(segments) %>%
        rename_(.dots=setNames("seqnames", "chromosome")) %>%
        select_(~chromosome, ~start, ~end)
    write.table(seg.df, file=file, quote=FALSE, sep="\t", row.names=FALSE)
}
