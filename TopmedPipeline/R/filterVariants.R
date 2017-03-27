#' Filter SeqVarGDSClass objects
#'
#' Set variant filters on SeqVarGDSClass objects
#'
#' These functions make it easy to apply various types of filters commonly used in association test code.
#' Most functions will not override previous filters, so they can be used in succession to apply multiple filters.
#' The one exception is \code{filterBySegment}, which must always be applied first.
#'
#' @param gds A \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @param segment An integer indicating which segment to select
#' @param segment.file The name of the file describing segments (see \code{\link{writeSegmentFile}})
#' @param pad.right The number of bases to add to the right of the segment (useful for sliding windows)
#' @param verbose Logical for whether to print number of variants selected, etc.
#' @name filterVariants
filterBySegment <- function(gds, segment, segment.file, pad.right=0, verbose=TRUE) {
    segments <- getSegments(segment.file)
    seg <- segments[segment]
    if (pad.right > 0) {
        seg <- resize(seg, width(seg) + pad.right, fix="start")
    }
    seqSetFilter(gds, variant.sel=seg, verbose=verbose)
}

#' @param idfile RData file with vector of variant.id
#' @rdname filterVariants
filterByFile <- function(gds, idfile, verbose=TRUE) {
    variant.id <- getobj(idfile)
    seqSetFilter(gds, variant.id=variant.id, action="intersect", verbose=verbose)
}

#' @param chr Chromosome to select
#' @rdname filterVariants
filterByChrom <- function(gds, chr, verbose=TRUE) {
    chrom <- seqGetData(gds, "chromosome")
    seqSetFilter(gds, variant.sel=(chrom == chr), action="intersect", verbose=verbose)
}

#' @rdname filterVariants
filterByPass <- function(gds, verbose=TRUE) {
    filt <- seqGetData(gds, "annotation/filter")
    seqSetFilter(gds, variant.sel=(filt == "PASS"), action="intersect", verbose=verbose)
}

#' @param biallelic Logical for whether to select only biallelic SNVs
#' @rdname filterVariants
filterBySNV <- function(gds, biallelic=TRUE, verbose=TRUE) {
    snv <- isSNV(gds, biallelic=biallelic)
    seqSetFilter(gds, variant.sel=snv, action="intersect", verbose=verbose)
}

#' @param sample.id Samples to include in calculating MAF
#' @param mac.min Minimum minor allele count (effective N) to include
#' @param maf.min Minimum MAF to include
#' @rdname filterVariants
filterByMAF <- function(gds, sample.id=NULL, mac.min=NA, maf.min=NA, verbose=TRUE) {
    if (sum(seqGetFilter(gds)$variant.sel) == 0) return(invisible())
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

#' @param build Genome build to use when identifying regions to exclude from PCA because of high correlation (HLA, LCT, inversions)
#' @rdname filterVariants
filterByPCAcorr <- function(gds, build="hg19", verbose=TRUE) {
    filt <- get(data(list=paste("pcaSnpFilters", build, sep="."), package="GWASTools"))
    chrom <- seqGetData(gds, "chromosome")
    pos <- seqGetData(gds, "position")
    pca.filt <- rep(TRUE, length(chrom))
    for (f in 1:nrow(filt)) {
        pca.filt[chrom == filt$chrom[f] & filt$start.base[f] < pos & pos < filt$end.base[f]] <- FALSE
    }
    seqSetFilter(gds, variant.sel=pca.filt, action="intersect", verbose=verbose)
}


#' @rdname filterVariants
checkSelectedVariants <- function(gds) {
    nvar <- sum(seqGetFilter(gds)$variant.sel)
    if (nvar == 0) {
        message("No variants selected. Exiting gracefully.")
        q(save="no", status=0)
    } else {
        message("Selected ", nvar, " variants.")
    }
}
