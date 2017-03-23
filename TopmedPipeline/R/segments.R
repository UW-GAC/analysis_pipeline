defineSegments <- function(seg.length=1e7, build="hg19") {
    gr <- get(data(list=paste("chromosomes", build, sep="_"),
                   package="TopmedPipeline", envir=environment()))
    do.call(c, lapply(gr, function(x) {
        window.start <- seq(BiocGenerics::start(x), BiocGenerics::end(x), seg.length)
        window.end <- seq(BiocGenerics::start(x) + seg.length - 1, BiocGenerics::end(x) + seg.length, by=seg.length)
        GRanges(seqnames=seqnames(x), IRanges(window.start, window.end))
    }))
}

writeSegmentFile <- function(segments, file) {
    seg.df <- as.data.frame(segments) %>%
        rename_(.dots=setNames("seqnames", "chromosome")) %>%
        select_(~chromosome, ~start, ~end)
    write.table(seg.df, file=file, quote=FALSE, sep="\t", row.names=FALSE)
}

getSegments <- function(file) {
    dat <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    GRanges(seqnames=dat$chromosome,
            ranges=IRanges(start=dat$start, end=dat$end))
}

# return the elements of varList where the first variant is in the segment
subsetBySegment <- function(varList, segment, segment.file) {
    # create a GRanges object containing the first variant from each item in varList
    dat <- do.call(rbind, lapply(varList, function(x) x[1,]))
    gr <- GRanges(seqnames=dat$chromosome,
                  ranges=IRanges(start=dat$position, end=dat$position))
    
    segments <- getSegments(segment.file)
    ind <- queryHits(findOverlaps(gr, segments[segment]))
    varList[ind]
}
