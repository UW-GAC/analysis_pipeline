filterByChrom <- function(gds, chr) {
    variant.id <- seqGetData(gds, "variant.id")
    chrom <- seqGetData(gds, "chromosome")
    variant.id <- variant.id[chrom == chr]
    seqSetFilter(gds, variant.id=variant.id)
    gds
}

filterByPass <- function(gds) {
    variant.id <- seqGetData(gds, "variant.id")
    filt <- seqGetData(gds, "annotation/filter")
    variant.id <- variant.id[filt == "PASS"]
    seqSetFilter(gds, variant.id=variant.id)
    gds
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
        variant.id <- seqGetData(gds, "variant.id")
        variant.id <- variant.id[maf.filt]
        seqSetFilter(gds, variant.id=variant.id)
    }
    gds
}
