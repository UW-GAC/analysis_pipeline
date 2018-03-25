context("variantData tests")
library(SeqVarTools)

.opengds <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    SeqVarData(seqOpen(seqExampleFileName()))
}

.testVariantsID <- function(gds, n=100) {
    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=sort(sample(1:SeqVarTools:::.nVar(gds), n)), verbose=FALSE)
    dat <- data.frame(variant.id=seqGetData(gds, "variant.id"),
                      weight=rnorm(100, mean=1, sd=0.1),
                      stringsAsFactors=FALSE)
    seqResetFilter(gds, verbose=FALSE)

    dat
}

.testVariantsPos <- function(gds, n=100) {
    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=sort(sample(1:SeqVarTools:::.nVar(gds), n)), verbose=FALSE)
    dat <- data.frame(chr=seqGetData(gds, "chromosome"),
                      pos=seqGetData(gds, "position"),
                      ref=refChar(gds),
                      alt=altChar(gds),
                      weight=rnorm(100, mean=1, sd=0.1),
                      stringsAsFactors=FALSE)
    seqResetFilter(gds, verbose=FALSE)

    dat
}

test_that("error for wrong columns in data.frame", {
    gds <- .opengds()
    expect_error(addVariantData(gds, variants=data.frame(a=100)))
    seqClose(gds)
})

test_that("match on variant.id", {
    gds <- .opengds()
    n <- 100
    var <- .testVariantsID(gds, n=n)
    gds <- addVariantData(gds, variants=var)
    weight <- variantData(gds)$weight
    expect_equal(sum(!is.na(weight)), n)
    expect_equivalent(na.omit(weight), var$weight)
    seqClose(gds)
})

test_that("match on pos", {
    gds <- .opengds()
    n <- 100
    var <- .testVariantsPos(gds, n=n)
    gds <- addVariantData(gds, variants=var)
    weight <- variantData(gds)$weight
    expect_equal(sum(!is.na(weight)), n)
    expect_equivalent(na.omit(weight), var$weight)
    seqClose(gds)
})
