context("utils tests")

test_that("constructFilename", {
    expect_equal(constructFilename("a", chromosome="1", segment="2"), "a_chr1_seg2.RData")
    expect_equal(constructFilename("a", chromosome="1"), "a_chr1.RData")
    expect_equal(constructFilename("a", segment="2"), "a_seg2.RData")
    expect_equal(constructFilename("a"), "a.RData")
})
