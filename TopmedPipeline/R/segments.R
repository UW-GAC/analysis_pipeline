#' Start and end positions of each chromosome
#'
#' @format \code{\link[GenomicRanges]{GRanges}} object
#' @source UCSC Genome Browser
#' @examples
#' \dontrun{
#' library(dplyr)
#' db <- src_mysql(dbname="hg19", host="genome-mysql.cse.ucsc.edu", user="genome")
#' chromInfo <- tbl(db, "chromInfo")
#' 
#' dat <- select(chromInfo, chrom, size) %>%
#'     collect() %>%
#'     filter(chrom %in% paste0("chr", c(1:22, "X"))) %>%
#'     mutate(chrom=factor(sub("chr", "", chrom), levels=c(1:22, "X"))) %>%
#'     arrange(chrom) %>%
#'     mutate(chrom=as.character(chrom))
#' 
#' library(GenomicRanges)
#' chromosomes_hg19 <- GRanges(seqnames=dat$chrom,
#'                             ranges=IRanges(start=1, end=dat$size))
#' }
"chromosomes_hg19"

#' @rdname chromosomes_hg19
"chromosomes_hg38"


#' Define segments for parallel processing of the genome
#'
#' @param seg.length Segment length in base pairs
#' @param n Number of segments to generate. Only used if \code{seg.length} is missing.
#'   Minimum number of segments is 23 (one per chromosome).
#' @param build Genome build
#' @return \code{\link[GenomicRanges]{GRanges}} object with segments covering the genome
#' @seealso \code{\link{writeSegmentFile}}, \code{\link{subsetBySegment}}
#' 
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom IRanges IRanges
#' @importFrom utils data
#' @export
defineSegments <- function(seg.length, n, build=c("hg19", "hg38")) {
    # load GRanges object with chromosomes for this build
    build <- match.arg(build)
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
    do.call(c, lapply(seq_along(gr), function(i) {
        x <- gr[i]
        window.start <- seq(BiocGenerics::start(x), BiocGenerics::end(x), x$seg.length)
        GRanges(seqnames=seqnames(x), IRanges(window.start, width=x$seg.length))
    }))
}

#' Write segments to a text file
#'
#' @param segments \code{\link[GenomicRanges]{GRanges}} object with segments
#' @param file Output file name
#' @seealso \code{\link{defineSegments}}, \code{\link{getSegments}}, \code{\link{filterBySegment}}
#'
#' @importFrom dplyr "%>%" rename_ select_
#' @importFrom stats setNames
#' @importFrom utils write.table
#' @export
writeSegmentFile <- function(segments, file) {
    seg.df <- as.data.frame(segments) %>%
        rename_(.dots=setNames("seqnames", "chromosome")) %>%
        select_(~chromosome, ~start, ~end)
    write.table(seg.df, file=file, quote=FALSE, sep="\t", row.names=FALSE)
}

#' Read segments from a text file
#'
#' @param file File with column names "chromosome", "start", "end " in the header
#' @return \code{\link[GenomicRanges]{GRanges}} object with segments
#' @seealso \code{\link{defineSegments}}, \code{\link{writeSegmentFile}}, \code{\link{subsetBySegment}}
#' 
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom utils read.table
#' @export
getSegments <- function(file) {
    dat <- read.table(file, header=TRUE, stringsAsFactors=FALSE)
    stopifnot(setequal(names(dat), c("chromosome", "start", "end")))
    GRanges(seqnames=dat$chromosome,
            ranges=IRanges(start=dat$start, end=dat$end))
}


#' Subset list of variant groups by segment
#' 
#' Select groups of variants where the first variant in the group is in the requested segment
#'
#' @param varList A list of data.frames, each with columns "chr" and "pos" OR a GRanges or GRangesList
#' @param segment An integer indicating which segment to select
#' @param segment.file The name of the file describing segments
#' @return Subset of \code{varList} where the first variant is in the segment \code{segment}
#' @seealso \code{\link{defineSegments}}, \code{\link{writeSegmentFile}}, \code{\link{filterBySegment}}, \code{\link{aggregateGRanges}}
#' 
#' @importFrom GenomicRanges findOverlaps GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors endoapply queryHits
#' @importFrom methods is
#' @export
subsetBySegment <- function(varList, segment, segment.file) {
    # create a GRanges object containing the first variant from each item in varList
    if (is(varList, "GRangesList")) {
        gr <- BiocGenerics::unlist(endoapply(varList, function(x) x[1]))
    } else if (is(varList, "GRanges")) {
        gr <- varList
    } else if (is.list(varList)) {
        dat <- do.call(rbind, lapply(varList, function(x) x[1,]))
        gr <- GRanges(seqnames=dat$chr,
                      ranges=IRanges(start=dat$pos, end=dat$pos))
    } else {
        stop("varList must be a list, GRanges, or GRangesList")
    }
    
    segments <- getSegments(segment.file)
    ind <- queryHits(findOverlaps(gr, segments[segment]))
    varList[ind]
}
