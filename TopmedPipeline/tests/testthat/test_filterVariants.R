context("filterVariants tests")
library(gdsfmt)
library(GenomicRanges)

.testData <- function() {
    showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- seqExampleFileName("gds")
    gds <- seqOpen(gdsfile)
}


test_that("filterBySegment", {
    gds <- .testData()
    data(segments)
    gr <- granges(gds)
    ol <- findOverlaps(gr, segments[1])
    filterBySegment(gds, 1, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), length(queryHits(ol)))

    seqClose(gds)
})


test_that("filterByFile", {
    gds <- .testData()
    id <- seqGetData(gds, "variant.id")[1:100]
    idfile <- tempfile()
    save(id, file=idfile)
    filterByFile(gds, idfile, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), 100)

    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=1:10, verbose=FALSE)
    filterByFile(gds, idfile, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), 10)

    seqClose(gds)
    unlink(idfile)
})


test_that("filterByChrom", {
    gds <- .testData()
    chr <- seqGetData(gds, "chromosome")
    filterByChrom(gds, 1, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), sum(chr == 1))

    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=1:10, verbose=FALSE)
    filterByChrom(gds, 1, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), 10)

    seqClose(gds)
})


test_that("filterByPass", {
    gds <- .testData()
    filt <- seqGetData(gds, "annotation/filter")
    filterByPass(gds, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), sum(filt == "PASS"))

    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=1:10, verbose=FALSE)
    filterByPass(gds, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), 10)

    seqClose(gds)
})


test_that("filterByMAF", {
    gds <- .testData()
    freq <- seqAlleleFreq(gds)
    maf <- pmin(freq, 1-freq)
    filterByMAF(gds, mac.min=NA, maf.min=0.1, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), sum(maf >= 0.1))

    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=1:10, verbose=FALSE)
    filterByMAF(gds, mac.min=NA, maf.min=0.1, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), sum(maf[1:10] >= 0.1))

    seqClose(gds)
})


