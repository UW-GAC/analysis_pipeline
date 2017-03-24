defineSegments <- function(seg.length, n, build="hg19") {
    # load GRanges object with chromosomes for this build
    gr <- get(data(list=paste("chromosomes", build, sep="_"),
                   package="TopmedPipeline", envir=environment()))

    # define segment length
    if (missing(seg.length)) {
        if (missing(n)) stop("must supply either seg.length or n")
        genome.length <- sum(as.numeric(width(gr)))
        genome.frac <- as.numeric(width(gr)) / genome.length
        n.chr <- round(genome.frac * n)
        if (sum(n.chr) > n) n.chr <- floor(genome.frac * n)
        if (any(n.chr == 0)) n.chr <- rep(1, length(gr))
        gr$seg.length <- ceiling(width(gr)/n.chr)
    } else {
        gr$seg.length <- seg.length
    }

    # create segments
    do.call(c, lapply(gr, function(x) {
        window.start <- seq(BiocGenerics::start(x), BiocGenerics::end(x), x$seg.length)
        GRanges(seqnames=seqnames(x), IRanges(window.start, width=x$seg.length))
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
