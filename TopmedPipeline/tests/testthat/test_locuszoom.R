context("LocusZoom tests")

test_that("calculateLD", {
    gdsfile <- seqExampleFileName("gds")
    r <- suppressWarnings(calculateLD(gdsfile, variant.id=1:10))
    expect_equal(dim(r), c(10,10))
    r <- suppressWarnings(calculateLD(gdsfile, variant.id=1:10, ref.var=10))
    expect_equal(length(r), 10)
})
