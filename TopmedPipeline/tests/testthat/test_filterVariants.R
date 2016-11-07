context("filterVariants tests")
library(gdsfmt)

.testData <- function() {
    showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- seqExampleFileName("gds")
    gds <- seqOpen(gdsfile)
}

test_that("filterByChrom", {
    gds <- .testData()
    chr <- seqGetData(gds, "chromosome")
    filterByChrom(gds, 1)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), sum(chr == 1))

    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=1:10)
    filterByChrom(gds, 1)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), 10)

    seqClose(gds)
})


test_that("filterByPass", {
    gds <- .testData()
    filt <- seqGetData(gds, "annotation/filter")
    filterByPass(gds)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), sum(filt == "PASS"))

    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=1:10)
    filterByPass(gds)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), 10)

    seqClose(gds)
})


test_that("filterByMAF", {
    gds <- .testData()
    freq <- seqAlleleFreq(gds)
    maf <- pmin(freq, 1-freq)
    sample.id <- seqGetData(gds, "sample.id")
    filterByMAF(gds, sample.id=sample.id, mac.min=NA, maf.min=0.1)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), sum(maf >= 0.1))

    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=1:10)
    filterByMAF(gds, sample.id=sample.id, mac.min=NA, maf.min=0.1)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), sum(maf[1:10] >= 0.1))

    seqClose(gds)
})


