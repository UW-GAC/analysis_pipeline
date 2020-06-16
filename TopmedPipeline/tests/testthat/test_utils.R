context("utils tests")

test_that("constructFilename", {
    expect_equal(constructFilename("a", chromosome="1", segment="2"), "a_chr1_seg2.RData")
    expect_equal(constructFilename("a", chromosome="1"), "a_chr1.RData")
    expect_equal(constructFilename("a", segment="2"), "a_seg2.RData")
    expect_equal(constructFilename("a"), "a.RData")
})

test_that("list2gds", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    x <- list(a=letters,
              b=1:10,
              c=rnorm(10),
              d=matrix(sample(1:100), nrow=10, ncol=10))
    f <- tempfile()
    list2gds(x, f)
    gds <- gdsfmt::openfn.gds(f)
    for (v in names(x)) {
        expect_equal(x[[v]], gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, v)))
    }
    gdsfmt::closefn.gds(gds)
    unlink(f)
})

test_that("gds2ibdobj", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gds <- seqOpen(seqExampleFileName())
    ibd <- SNPRelate::snpgdsIBDKING(gds, verbose=FALSE)
    gdsfile <- tempfile()
    list2gds(ibd, gdsfile)
    ibd$afreq <- NULL
    ibd2 <- gds2ibdobj(gdsfile)
    expect_equal(ibd, ibd2)

    samp.sel <- sort(sample(seq_along(ibd$sample.id), 10))
    sel <- SNPRelate::snpgdsIBDSelection(ibd, samp.sel=samp.sel)
    ibd2.sel <- gds2ibdobj(gdsfile, sample.id=ibd2$sample.id[samp.sel])
    sel2 <- SNPRelate::snpgdsIBDSelection(ibd2.sel)
    expect_equal(sel, sel2)
})

test_that("rds", {
    rdsfile <- paste0(tempfile(), ".rds")
    x <- 1:10
    saveRDS(x, rdsfile)
    x2 <- getobj(rdsfile)
    expect_equal(x, x2)
    unlink(rdsfile)
})
