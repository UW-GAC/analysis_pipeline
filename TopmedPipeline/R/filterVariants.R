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
    #ind
}

## always set this one first (no intersect option)
filterBySegment <- function(gds, segment, segment.file, verbose=TRUE) {
    segments <- getSegments(segment.file)
    seqSetFilter(gds, variant.sel=segments[segment], verbose=verbose)
}

filterByFile <- function(gds, idfile, verbose=TRUE) {
    variant.id <- getobj(idfile)
    seqSetFilter(gds, variant.id=variant.id, action="intersect", verbose=verbose)
}

filterByChrom <- function(gds, chr, verbose=TRUE) {
    chrom <- seqGetData(gds, "chromosome")
    seqSetFilter(gds, variant.sel=(chrom == chr), action="intersect", verbose=verbose)
}

filterByPass <- function(gds, verbose=TRUE) {
    filt <- seqGetData(gds, "annotation/filter")
    seqSetFilter(gds, variant.sel=(filt == "PASS"), action="intersect", verbose=verbose)
}

filterByMAF <- function(gds, sample.id=NULL, mac.min=NA, maf.min=NA, verbose=TRUE) {
    if ((!is.na(mac.min) & mac.min > 1) |
        (!is.na(maf.min) & maf.min > 0)) {
        if (is.null(sample.id)) sample.id <- seqGetData(gds, "sample.id")
        seqSetFilter(gds, sample.id=sample.id, verbose=FALSE)
        ref.freq <- seqAlleleFreq(gds)
        maf <- pmin(ref.freq, 1-ref.freq)
        if (!is.na(mac.min)) {
            maf.filt <- 2 * maf * (1-maf) * length(sample.id) >= mac.min
            if (verbose) message(paste("Running on", sum(maf.filt), "variants with MAC >=", mac.min))
        } else {
            maf.filt <- maf >= maf.min
            if (verbose) message(paste("Running on", sum(maf.filt), "variants with MAF >=", maf.min))
        }
        seqSetFilter(gds, variant.sel=maf.filt, action="intersect", verbose=verbose)
    }
}

checkSelectedVariants <- function(gds) {
    nvar <- sum(seqGetFilter(gds)$variant.sel)
    if (nvar == 0) {
        message("No variants selected. Exiting gracefully.")
        q(save="no", status=0)
    } else {
        message("Selected ", nvar, " variants.")
    }
}
