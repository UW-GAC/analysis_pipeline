context("thinPoints tests")

test_that("no grouping", {
    dat <- data.frame(x=1:100)
    thin <- thinPoints(dat, "x", n=2, nbins=10, groupBy=NULL)
    expect_equal(nrow(thin), 20)
    expect_true(all(table(cut(thin$x, breaks=seq(0,100,10), labels=FALSE)) == 2))
})

test_that("with grouping", {
    dat <- data.frame(x=1:200, y=c(rep("a",100), rep("b",100)))
    thin <- thinPoints(dat, "x", n=2, nbins=10, groupBy="y")
    expect_equal(nrow(thin), 40)
    expect_true(all(table(cut(thin$x[thin$y == "a"], breaks=seq(0,100,10), labels=FALSE)) == 2))
    expect_true(all(table(cut(thin$x[thin$y == "b"], breaks=seq(0,100,10), labels=FALSE)) == 2))
})

test_that("bins with small n are kept", {
    dat <- data.frame(x=c(1:50, 99:100))
    thin <- thinPoints(dat, "x", n=10, nbins=2)
    expect_true(all(99:100 %in% thin$x))
})
