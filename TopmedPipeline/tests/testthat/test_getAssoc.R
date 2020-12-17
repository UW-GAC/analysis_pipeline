context("getAssoc tests")
library(dplyr)
library(GENESIS)
library(gdsfmt)
library(SeqVarTools)
library(Biobase)
library(GenomicRanges)

.testData <- function() {
    showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- seqExampleFileName("gds")
    gds <- seqOpen(gdsfile)

    data(sample_annotation)
    SeqVarData(gds, sampleData=sample_annotation)
}

.testNullModel <- function(seqData, MM=FALSE) {
    if (MM) {
        data(grm)
    } else {
        grm <- NULL
    }
    fitNullModel(sampleData(seqData), outcome="outcome", covars="sex", cov.mat=grm, verbose=FALSE)
}


test_that("single related", {
    seqData <- SeqVarBlockIterator(.testData(), verbose=FALSE)
    nullmod <- .testNullModel(seqData, MM=TRUE)
    seqSetFilterChrom(seqData, include=1, verbose=FALSE)
    a1 <- assocTestSingle(seqData, nullmod, verbose=FALSE)
    seqSetFilterChrom(seqData, include=2, verbose=FALSE)
    a2 <- assocTestSingle(seqData, nullmod, verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])
    a <- rbind(a1, a2) %>%
        filter_(~!is.na(Score.pval))

    assoc <- getAssoc(files, "single")
    expect_equal(nrow(assoc), nrow(a))
    expect_equal(assoc$id, a$variant.id)
    expect_equal(as.character(assoc$chr), a$chr)
    expect_equal(assoc$pos, a$pos)
    expect_equal(assoc$start, assoc$pos)
    expect_equal(assoc$end, assoc$pos)
    expect_equal(assoc$stat, a$Score.Stat)
    expect_equal(assoc$pval, a$Score.pval)

    seqClose(seqData)
    unlink(files)
})


test_that("window, burden", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    seqSetFilterChrom(seqData, include=1, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a1 <- assocTestAggregate(seqData, nullmod, verbose=FALSE)
    seqSetFilterChrom(seqData, include=2, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a2 <- assocTestAggregate(seqData, nullmod, verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])
    a <- rbind(a1$results, a2$results) %>%
        filter(n.site > 0) %>%
        mutate(id = sprintf("chr%s_%d", chr, start))

    assoc <- getAssoc(files, "window")
    expect_equal(nrow(assoc), nrow(a))
    expect_equal(assoc$id, a$id)
    expect_equal(assoc$stat, a$Score.Stat)
    expect_equal(assoc$pval, a$Score.pval)

    seqClose(seqData)
    unlink(files)
})


test_that("aggregate, skat", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    gr <- granges(seqData)
    agg <- GRangesList(gr[1:100], gr[101:200], gr[201:300])
    names(agg) <- letters[1:3]
    seqData <- SeqVarListIterator(seqData, variantRanges=agg, verbose=FALSE)
    a1 <- assocTestAggregate(seqData, nullmod, test="SKAT", verbose=FALSE)
    seqResetFilter(seqData, verbose=FALSE)
    agg <- GRangesList(gr[301:400], gr[401:500], gr[501:600])
    names(agg) <- letters[4:6]
    seqData <- SeqVarListIterator(seqData, variantRanges=agg, verbose=FALSE)
    a2 <- assocTestAggregate(seqData, nullmod, test="SKAT", verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])
    a <- rbind(a1$results, a2$results) %>%
        filter(n.site > 0) %>%
        mutate(id = c(names(a1$variantInfo), names(a2$variantInfo)))

    assoc <- getAssoc(files, "aggregate")
    expect_equal(nrow(assoc), nrow(a))
    expect_equal(assoc$id, a$id)
    expect_true(setequal(assoc$pval, a$pval))
    expect_equal(as.character(assoc$chr[1]), a1$variantInfo[[1]]$chr[1])
    expect_true(assoc$pos[1] > a1$variantInfo[[1]]$pos[1] & assoc$pos[1] < max(a1$variantInfo[[1]]$pos))
    expect_equal(assoc$start[1], a1$variantInfo[[1]]$pos[1])
    expect_equal(assoc$end[1], max(a1$variantInfo[[1]]$pos))

    seqClose(seqData)
    unlink(files)
})

test_that("aggregate, multiple chromosomes in one aggregation unit", {
  seqData <- .testData()
  nullmod <- .testNullModel(seqData)
  gr <- granges(seqData)
  gr_chr1 <- gr[(seqnames(gr) == 1)][1:10]
  gr_chr2 <- gr[(seqnames(gr) == 2)]
  agg <- GRangesList(c(gr_chr1, gr_chr2))
  names(agg) <- "a"
  seqData <- SeqVarListIterator(seqData, variantRanges=agg, verbose=FALSE)
  a <- assocTestAggregate(seqData, nullmod, test="SKAT", verbose=FALSE)
  file <- tempfile()
  save(a, file = file)

  assoc <- getAssoc(file, "aggregate")

  tmp <- a$variantInfo[[1]] %>%
    filter(chr == 2)
  expected_pos <- floor((min(tmp$pos) + max(tmp$pos)) / 2)
  expect_equal(nrow(assoc), 1)
  expect_equal(assoc$id, "a")
  expect_equal(as.character(assoc$chr), "2")
  expect_equal(assoc$pos, expected_pos)

  # Chooses first chromosome if there is an equal number of variants on each chromosome.
  a$variantInfo[[1]] <- a$variantInfo[[1]] %>%
    group_by(chr) %>%
    sample_n(2) %>%
    ungroup()
  save(a, file = file)

  tmp <- a$variantInfo[[1]] %>%
    filter(chr == chr[1])
  expected_pos <- floor((min(tmp$pos) + max(tmp$pos)) / 2)
  assoc <- getAssoc(file, "aggregate")
  expect_equal(nrow(assoc), 1)
  expect_equal(assoc$id, "a")
  expect_equal(as.character(assoc$chr), "1")
  expect_equal(assoc$pos, expected_pos)

})

test_that("window, smmat", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    seqSetFilterChrom(seqData, include=1, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a1 <- assocTestAggregate(seqData, nullmod, test="SMMAT", verbose=FALSE)
    seqSetFilterChrom(seqData, include=2, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a2 <- assocTestAggregate(seqData, nullmod, test="SMMAT", verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])
    a <- rbind(a1$results, a2$results) %>%
        filter_(~(n.site > 0))

    assoc <- getAssoc(files, "window")
    expect_equal(assoc$pval, a$pval_SMMAT)

    seqClose(seqData)
    unlink(files)
})


test_that("window, skato", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    seqSetFilterChrom(seqData, include=1, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a1 <- assocTestAggregate(seqData, nullmod, test="SKATO", rho=c(0,1), verbose=FALSE)
    seqSetFilterChrom(seqData, include=2, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a2 <- assocTestAggregate(seqData, nullmod, test="SKATO", rho=c(0,1), verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])
    a <- rbind(a1$results, a2$results) %>%
        filter_(~(n.site > 0))

    assoc <- getAssoc(files, "window")
    expect_equal(assoc$pval, a$pval_SKATO)

    seqClose(seqData)
    unlink(files)
})


test_that("combine single", {
    seqData <- SeqVarBlockIterator(.testData(), verbose=FALSE)
    nullmod <- .testNullModel(seqData, MM=TRUE)
    seqSetFilterChrom(seqData, include=1, verbose=FALSE)
    a1 <- assocTestSingle(seqData, nullmod, verbose=FALSE)
    seqSetFilterChrom(seqData, include=2, verbose=FALSE)
    a2 <- assocTestSingle(seqData, nullmod, verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])

    assoc <- combineAssoc(files, "single")

    seqSetFilterChrom(seqData, include=1:2, verbose=FALSE)
    a <- assocTestSingle(seqData, nullmod, verbose=FALSE)
    expect_true(all(a == assoc, na.rm=TRUE))
    expect_true(all(is.na(a) == is.na(assoc)))

    seqClose(seqData)
    unlink(files)
})

test_that("combine window", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    seqSetFilterChrom(seqData, include=1, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a1 <- assocTestAggregate(seqData, nullmod, verbose=FALSE)
    seqSetFilterChrom(seqData, include=2, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a2 <- assocTestAggregate(seqData, nullmod, verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])

    assoc <- combineAssoc(files, "window")

    seqData <- .testData()
    seqSetFilterChrom(seqData, include=1:2, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a <- assocTestAggregate(seqData, nullmod, verbose=FALSE)
    expect_equal(a, assoc)

    seqClose(seqData)
    unlink(files)
})

test_that("combine aggregate", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    gr <- granges(seqData)
    agg1 <- GRangesList(gr[1:100], gr[101:200], gr[201:300])
    seqData <- SeqVarListIterator(seqData, variantRanges=agg1, verbose=FALSE)
    a1 <- assocTestAggregate(seqData, nullmod, test="SKAT", verbose=FALSE)
    seqResetFilter(seqData, verbose=FALSE)
    agg2 <- GRangesList(gr[301:400], gr[401:500], gr[501:600])
    seqData <- SeqVarListIterator(seqData, variantRanges=agg2, verbose=FALSE)
    a2 <- assocTestAggregate(seqData, nullmod, test="SKAT", verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])

    assoc <- combineAssoc(files, "aggregate")

    seqResetFilter(seqData, verbose=FALSE)
    seqData <- SeqVarListIterator(seqData, variantRanges=c(agg1, agg2), verbose=FALSE)
    a <- assocTestAggregate(seqData, nullmod, test="SKAT", verbose=FALSE)
    expect_equal(a, assoc)

    seqClose(seqData)
    unlink(files)
})

test_that("combine aggregate - empty set", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    gr <- granges(seqData)
    agg1 <- GRangesList(GRanges(seqnames=1, IRanges(start=2000000, width=1)))
    seqData <- SeqVarListIterator(seqData, variantRanges=agg1, verbose=FALSE)
    a1 <- assocTestAggregate(seqData, nullmod, test="SKAT", verbose=FALSE)
    seqResetFilter(seqData, verbose=FALSE)
    agg2 <- GRangesList(gr[1:100], gr[101:200], gr[201:300])
    seqData <- SeqVarListIterator(seqData, variantRanges=agg2, verbose=FALSE)
    a2 <- assocTestAggregate(seqData, nullmod, test="SKAT", verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])

    assoc <- combineAssoc(files, "aggregate")

    seqResetFilter(seqData, verbose=FALSE)
    seqData <- SeqVarListIterator(seqData, variantRanges=c(agg1, agg2), verbose=FALSE)
    a <- assocTestAggregate(seqData, nullmod, test="SKAT", verbose=FALSE)
    expect_equal(a, assoc)

    seqClose(seqData)
    unlink(files)
})

test_that("arrange_chr_pos", {
    x <- data.frame(chr=c(10,2,"X",1,1), pos=c(1,1,1,2,1), stringsAsFactors=FALSE)
    xa <- x[c(5,4,2,1,3),]; rownames(xa) <- 1:5
    expect_equal(.arrange_chr_pos(x), xa)
    names(x) <- names(xa) <- c("a", "b")
    expect_equal(.arrange_chr_pos(x, chr="a", pos="b"), xa)
})

test_that("index_chr_pos", {
    x <- list(data.frame(chr=2, pos=1:5),
              data.frame(chr=1, pos=6:10),
              data.frame(chr="X", pos=1:5, stringsAsFactors=FALSE),
              data.frame(chr=1, pos=1:5))
    expect_equal(.index_chr_pos(x), c(4,2,1,3))
    xa <- lapply(x, function(xx) {names(xx) <- c("a", "b"); xx})
    expect_equal(.index_chr_pos(xa, chr="a", pos="b"), c(4,2,1,3))
})

test_that("combine out of order", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData, MM=TRUE)

    seqSetFilter(seqData, variant.sel=1:10, verbose=FALSE)
    seqData <- SeqVarBlockIterator(seqData, verbose=FALSE)
    a1 <- assocTestSingle(seqData, nullmod, verbose=FALSE)
    seqSetFilter(seqData, variant.sel=11:20, verbose=FALSE)
    seqData <- SeqVarBlockIterator(seqData, verbose=FALSE)
    a2 <- assocTestSingle(seqData, nullmod, verbose=FALSE)
    files <- file.path(tempdir(), c("a", "b"))
    save(a1, file=files[1])
    save(a2, file=files[2])
    assoc <- combineAssoc(files, "single", ordered=TRUE)
    a <- combineAssoc(files[c(2,1)], "single", ordered=TRUE)
    expect_equal(a, assoc)

    seqSetFilterChrom(seqData, include=1, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a1 <- assocTestAggregate(seqData, nullmod, verbose=FALSE)
    seqSetFilterChrom(seqData, include=2, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a2 <- assocTestAggregate(seqData, nullmod, verbose=FALSE)
    save(a1, file=files[1])
    save(a2, file=files[2])
    assoc <- combineAssoc(files, "window", ordered=TRUE)
    a <- combineAssoc(files[c(2,1)], "window", ordered=TRUE)
    expect_equal(a, assoc)

    seqResetFilter(seqData, verbose=FALSE)
    gr <- granges(seqData)
    agg1 <- GRangesList(gr[101:200], gr[1:100], gr[501:600])
    seqData <- SeqVarListIterator(seqData, variantRanges=agg1, verbose=FALSE)
    a1 <- assocTestAggregate(seqData, nullmod, test="SKAT", verbose=FALSE)
    seqResetFilter(seqData, verbose=FALSE)
    agg2 <- GRangesList(gr[301:400], gr[401:500], gr[201:300])
    seqData <- SeqVarListIterator(seqData, variantRanges=agg2, verbose=FALSE)
    a2 <- assocTestAggregate(seqData, nullmod, test="SKAT", verbose=FALSE)
    save(a1, file=files[1])
    save(a2, file=files[2])
    assoc <- combineAssoc(files, "aggregate", ordered=TRUE)
    rownames(assoc$results) <- 1:nrow(assoc$results)
    agg <- c(agg1, agg2)[c(2,1,6,4,5,3)]
    seqResetFilter(seqData, verbose=FALSE)
    seqData <- SeqVarListIterator(seqData, variantRanges=agg, verbose=FALSE)
    a <- assocTestAggregate(seqData, nullmod, test="SKAT", verbose=FALSE)
    expect_equal(a, assoc)

    seqClose(seqData)
    unlink(files)
})


test_that("omitKnownHits", {
    assoc <- data.frame(chr=c(rep(1,100), rep(2,100)),
                        pos=c(sort(sample(1:10000, 100)), sort(sample(1:10000, 100))))
    hits <- assoc[sample(1:nrow(assoc), 20),]
    res <- omitKnownHits(assoc, hits, flank=0)
    .pos <- function(x) paste0(x$chr, x$pos)
    expect_equal(setdiff(.pos(assoc), .pos(hits)), .pos(res))

    res <- omitKnownHits(assoc, hits, flank=1)
    expect_true(nrow(res) < nrow(assoc))
})



test_that("MAC - single", {
    seqData <- SeqVarBlockIterator(.testData(), verbose=FALSE)
    nullmod <- .testNullModel(seqData)
    seqSetFilterChrom(seqData, include=1, verbose=FALSE)
    a <- assocTestSingle(seqData, nullmod, verbose=FALSE)
    a <- addMAC(a, "single")
    expect_true("MAC" %in% names(a))

    files <- tempfile()
    save(a, file=files)
    assoc <- getAssoc(files, "single")
    expect_true("MAC" %in% names(assoc))

    seqClose(seqData)
    unlink(files)
})



test_that("MAC - window", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    seqSetFilterChrom(seqData, include=1, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a <- assocTestAggregate(seqData, nullmod, verbose=FALSE)
    a <- addMAC(a, "window")
    expect_true("MAC" %in% names(a$results))

    # for variants where alt is minor allele, n.alt should be MAC
    alt.min <- which(sapply(a$variantInfo, function(x) all(x$freq < 0.5)))
    expect_equal(a$results$n.alt[alt.min], a$results$MAC[alt.min])

    files <- tempfile()
    save(a, file=files)
    assoc <- getAssoc(files, "window")
    expect_true("MAC" %in% names(assoc))

    seqClose(seqData)
    unlink(files)
})


test_that("remove conditional", {
    seqData <- SeqVarBlockIterator(.testData(), verbose=FALSE)
    nullmod <- .testNullModel(seqData, MM=TRUE)
    seqSetFilterChrom(seqData, include=1, verbose=FALSE)
    a1 <- assocTestSingle(seqData, nullmod, verbose=FALSE)
    dat <- data.frame(variant.id=1:3,
                      chromosome=1,
                      stringsAsFactors=FALSE)
    cvfile <- tempfile()
    save(dat, file=cvfile)
    a2 <- removeConditional(a1, cvfile)
    expect_equivalent(a1[4:nrow(a1),], a2)
    unlink(cvfile)
})
