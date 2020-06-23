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
#'
#' @import SeqArray
#' @importFrom GenomicRanges resize width
#' @export
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
#'
#' @import SeqArray
#' @export
filterByFile <- function(gds, idfile, verbose=TRUE) {
    variant.id <- getobj(idfile)
    seqSetFilter(gds, variant.id=variant.id, action="intersect", verbose=verbose)
}

#' @param chr Chromosome to select
#' @rdname filterVariants
#'
#' @import SeqArray
#' @export
filterByChrom <- function(gds, chr, verbose=TRUE) {
    chrom <- seqGetData(gds, "chromosome")
    seqSetFilter(gds, variant.sel=(chrom == chr), action="intersect", verbose=verbose)
}

#' @rdname filterVariants
#'
#' @import SeqArray
#' @export
filterByPass <- function(gds, verbose=TRUE) {
    filt <- seqGetData(gds, "annotation/filter")
    seqSetFilter(gds, variant.sel=(filt == "PASS"), action="intersect", verbose=verbose)
}

#' @param biallelic Logical for whether to select only biallelic SNVs
#' @rdname filterVariants
#'
#' @import SeqArray
#' @importFrom SeqVarTools isSNV
#' @export
filterBySNV <- function(gds, biallelic=TRUE, verbose=TRUE) {
    snv <- isSNV(gds, biallelic=biallelic)
    seqSetFilter(gds, variant.sel=snv, action="intersect", verbose=verbose)
}

.addBuild <- function(gds, calc.fn, build, ...) {
    if (is(gds, "SeqVarData")) {
        return(calc.fn(gds, genome.build=build, ...))
    } else {
        return(calc.fn(gds, ...))
    }
}

#' @importFrom SeqVarTools alleleFrequency
.calcMAF <- function(gds, sample.id, build="hg19") {
    seqSetFilter(gds, sample.id=sample.id, verbose=FALSE)
    ref.freq <- .addBuild(gds, alleleFrequency, build)
    pmin(ref.freq, 1-ref.freq)
}

#' @importFrom SeqVarTools minorAlleleCount
.calcMAC <- function(gds, sample.id, build="hg19") {
    seqSetFilter(gds, sample.id=sample.id, verbose=FALSE)
    round(.addBuild(gds, minorAlleleCount, build))
}

#' @importFrom SeqVarTools alleleFrequency
.calcEffN <- function(gds, sample.id, build="hg19") {
    seqSetFilter(gds, sample.id=sample.id, verbose=FALSE)
    ref.freq <- .addBuild(gds, alleleFrequency, build)
    n.obs <- SeqVarTools:::.nSampObserved(gds)
    2 * n.obs * ref.freq * (1-ref.freq)
}

#' @importFrom SeqVarTools alleleFrequency
.calcAltFreq <- function(gds, sample.id, build="hg19") {
    seqSetFilter(gds, sample.id=sample.id, verbose=FALSE)
    .addBuild(gds, alleleFrequency, n=1, build=build)
}

#' @importFrom SeqVarTools sampleData
.calcBinary <- function(gds, sample.id, binary.outcome, calc.fn, ...) {
    outcome <- sampleData(gds)[[binary.outcome]]
    cases <- sample.id[outcome == 1]
    controls <- sample.id[outcome == 0]
    maf.list <- lapply(list(cases, controls), function(x) calc.fn(gds, sample.id=x, ...))
    pmin(maf.list[[1]], maf.list[[2]])
}

#' @param sample.id Samples to include in calculating allele frequency
#' @param maf.min Minimum MAF to include
#' @rdname filterVariants
#'
#' @import SeqArray
#' @export
filterByMAF <- function(gds, sample.id=NULL, maf.min=0, build="hg19", verbose=TRUE) {
    stopifnot(maf.min >= 0 & maf.min <= 0.5)
    if (maf.min == 0) return(invisible())
    if (sum(seqGetFilter(gds)$variant.sel) == 0) return(invisible())
    
    if (is.null(sample.id)) sample.id <- seqGetData(gds, "sample.id")
    maf <- .calcMAF(gds, sample.id, build=build)
    maf.filt <- maf >= maf.min
    if (verbose) message(paste("Running on", sum(maf.filt), "variants with MAF >=", maf.min))
    seqSetFilter(gds, variant.sel=maf.filt, action="intersect", verbose=verbose)
}

#' @param mac.min Minimum minor allele count to include
#' @rdname filterVariants
#' @export
filterByMAC <- function(gds, sample.id=NULL, mac.min=1, build="hg19", verbose=TRUE) {
    stopifnot(mac.min >= 0)
    if (mac.min == 0) return(invisible())
    if (sum(seqGetFilter(gds)$variant.sel) == 0) return(invisible())
    
    if (is.null(sample.id)) sample.id <- seqGetData(gds, "sample.id")
    mac <- .calcMAC(gds, sample.id, build=build)
    maf.filt <- mac >= mac.min
    if (verbose) message(paste("Running on", sum(maf.filt), "variants with MAC >=", mac.min))
    seqSetFilter(gds, variant.sel=maf.filt, action="intersect", verbose=verbose)
}

#' @param effN.min Minimum effective N, calculated as 2 * MAF * (1-MAF) * n.obs
#' @rdname filterVariants
#' @export
filterByEffN <- function(gds, sample.id=NULL, effN.min=1, build="hg19", verbose=TRUE) {
    stopifnot(effN.min >= 0)
    if (effN.min == 0) return(invisible())
    if (sum(seqGetFilter(gds)$variant.sel) == 0) return(invisible())
    
    if (is.null(sample.id)) sample.id <- seqGetData(gds, "sample.id")
    mac <- .calcEffN(gds, sample.id, build=build)
    maf.filt <- mac >= effN.min
    if (verbose) message(paste("Running on", sum(maf.filt), "variants with effN >=", effN.min))
    seqSetFilter(gds, variant.sel=maf.filt, action="intersect", verbose=verbose)
}

.minAltFreq <- function(f) {
    sapply(f, function(x) {
        x <- x[-1]
        x <- x[x > 0]
        if (length(x) > 0) return(min(x)) else return(NA)
    })
}

#' @param af.max Maximum alternate allele frequency to include
#' @rdname filterVariants
#'
#' @import SeqArray
#' @importFrom SeqVarTools alleleFrequency
#' @export
filterByRare <- function(gds, sample.id=NULL, af.max=0.01, build="hg19", verbose=TRUE) {
    stopifnot(af.max >= 0 & af.max <= 1)
    if (af.max == 1) return(invisible())
    if (sum(seqGetFilter(gds)$variant.sel) == 0) return(invisible())
    
    if (is.null(sample.id)) sample.id <- seqGetData(gds, "sample.id")
    alt.freq <- .calcAltFreq(gds, sample.id, build=build)
    af.filt <- alt.freq > 0 & alt.freq < 1 & alt.freq <= af.max

    ## check frequency for multiallelic variants
    n <- seqNumAllele(gds)
    multi <- which(n > 2)
    if (length(multi) > 0) {
        seqSetFilter(gds, variant.sel=multi, action="push+intersect", verbose=FALSE)
        # alleleFrequency requires selecting only one allele, so use seqAlleleFreq
        # this means we are skipping sex correction for multiallelic variants
        alt.freq <- seqAlleleFreq(gds, ref.allele=NULL)
        min.freq <- .minAltFreq(alt.freq)
        multi.filt <- !is.na(min.freq) & min.freq <= af.max
        af.filt[multi] <- af.filt[multi] | multi.filt
        seqSetFilter(gds, action="pop", verbose=FALSE)
    }
    
    if (verbose) message(paste("Running on", sum(af.filt), "variants with alternate allele frequency <=", af.max))
    seqSetFilter(gds, variant.sel=af.filt, action="intersect", verbose=verbose)
}

#' @param build Genome build to use when identifying regions to exclude from PCA because of high correlation (HLA, LCT, inversions)
#' @rdname filterVariants
#'
#' @import SeqArray
#' @importFrom utils data
#' @export
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
#'
#' @import SeqArray
#' @export
checkSelectedVariants <- function(gds) {
    nvar <- sum(seqGetFilter(gds)$variant.sel)
    if (nvar == 0) {
        message("No variants selected. Exiting gracefully.")
        q(save="no", status=0)
    } else {
        message("Selected ", nvar, " variants.")
    }
}
