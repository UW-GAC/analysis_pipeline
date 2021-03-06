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

test_that("calculateLambda", {
  # Null hypothesis.
  null_stat <- qchisq((0:1000) / 1000, 1)
  expect_equal(calculateLambda(null_stat, df = 1), 1)
  expect_equal(calculateLambda(null_stat * 0.9, df = 1), 0.9)
  expect_equal(calculateLambda(null_stat * 1.1, df = 1), 1.1)

  # Works with other quantiles
  expect_equal(calculateLambda(null_stat, df = 1, quantiles = 0.25), 1)
  expect_equal(calculateLambda(null_stat * 0.9, df = 1, quantiles = 0.25), 0.9)
  expect_equal(calculateLambda(null_stat * 1.1, df = 1, quantiles = 0.25), 1.1)
  expect_equal(calculateLambda(null_stat, df = 1, quantiles = 0.75), 1)
  expect_equal(calculateLambda(null_stat * 0.9, df = 1, quantiles = 0.75), 0.9)
  expect_equal(calculateLambda(null_stat * 1.1, df = 1, quantiles = 0.75), 1.1)

  # Works with multiple quantiles
  expect_equal(calculateLambda(null_stat, df = 1, quantiles = c(0.25, 0.5, 0.75)), c(1, 1, 1))
  expect_equal(calculateLambda(null_stat * 0.9, df = 1, quantiles = c(0.25, 0.5, 0.75)), c(0.9, 0.9, 0.9))
  expect_equal(calculateLambda(null_stat * 1.1, df = 1, quantiles = c(0.25, 0.5, 0.75)), c(1.1, 1.1, 1.1))
  # Check a test statistic that has different lambdas at different quantiles.
  idx_low <- 1:round(0.4 * length(null_stat))
  idx_high <- round(0.6 * length(null_stat)):length(null_stat)
  tmp_stat <- null_stat
  tmp_stat[idx_low] <- null_stat[idx_low] * 0.9
  tmp_stat[idx_high] <- null_stat[idx_high] * 1.1
  expect_equal(calculateLambda(tmp_stat, df = 1, quantiles = c(0.25, 0.5, 0.75)), c(0.9, 1, 1.1))

  # Largest stats are inflated.
  null_stat[800:999] <- null_stat[800:999]*2
  expect_equal(calculateLambda(null_stat, df = 1), 1)
  expect_equal(calculateLambda(null_stat, df = 1, quantiles = 0.05), 1)
  expect_equal(calculateLambda(null_stat, df = 1, quantiles = 0.95), 2)
  expect_equal(calculateLambda(null_stat, df = 1, quantiles = c(0.05, 0.5, 0.95)), c(1, 1, 2))

  # More degrees of freedom
  null_stat <- qchisq((0:1000) / 1000, 3)
  expect_equal(calculateLambda(null_stat, df = 3), 1)
  expect_equal(calculateLambda(null_stat * 0.9, df = 3), 0.9)
  expect_equal(calculateLambda(null_stat * 1.1, df = 3), 1.1)
})
