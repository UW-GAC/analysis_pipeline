filterByChrom <- function(gds, chr) {
    chrom <- seqGetData(gds, "chromosome")
    seqSetFilter(gds, variant.sel=(chrom == chr), action="intersect")
}

filterByPass <- function(gds) {
    filt <- seqGetData(gds, "annotation/filter")
    seqSetFilter(gds, variant.sel=(filt == "PASS"), action="intersect")
}

filterByMAF <- function(gds, sample.id, mac.min, maf.min) {
    if ((!is.na(mac.min) & mac.min > 1) |
        (!is.na(maf.min) & maf.min > 0)) {
        seqSetFilter(gds, sample.id=sample.id, verbose=FALSE)
        ref.freq <- seqAlleleFreq(gds)
        maf <- pmin(ref.freq, 1-ref.freq)
        if (!is.na(mac.min)) {
            maf.filt <- 2 * maf * (1-maf) * length(sample.id) >= mac.min
            message(paste("Running on", sum(maf.filt), "variants with MAC >=", mac.min))
        } else {
            maf.filt <- maf >= maf.min
            message(paste("Running on", sum(maf.filt), "variants with MAF >=", maf.min))
        }
        seqSetFilter(gds, variant.sel=maf.filt, action="intersect")
    }
}
