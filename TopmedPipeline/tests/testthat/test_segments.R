context("test splitting analyses by segment")
library(GENESIS)
library(gdsfmt)
library(SeqVarTools)
library(Biobase)
library(dplyr)
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

.testSegFile <- function(segments) {
    segfile <- tempfile()
    writeSegmentFile(segments, segfile)
    segfile
}

test_that("defineSegments", {
    data(chromosomes_hg19)
    sl <- sample(seq(1e6, 1e8, by=1000), 1)
    seg <- defineSegments(sl, build="hg19")
    expect_equal(reduce(seg) >= chromosomes_hg19, rep(TRUE, 23))
    seg2 <- defineSegments(sl, n=10, build="hg19") # n is ignored
    expect_equal(seg, seg2)

    n <- sample(100:1000, 1)
    seg <- defineSegments(n=n, build="hg19")
    expect_true(n - length(seg) < 23)
    expect_equal(reduce(seg) >= chromosomes_hg19, rep(TRUE, 23))

    n <- sample(1:23, 1)
    seg <- defineSegments(n=n, build="hg19")
    expect_true(length(seg) == 23)
    expect_equivalent(seg, chromosomes_hg19)
    
    expect_error(defineSegments())
})

test_that("read and write segment files", {
    data(segments)
    segfile <- tempfile()
    writeSegmentFile(segments, segfile)
    seg2 <- getSegments(segfile)
    expect_equivalent(segments, seg2)
    unlink(segfile)
})


test_that("single", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData, MM=TRUE)
    data(segments)
    segments <- segments[seqnames(segments) == 1]
    files <- character(length(segments))
    for (i in seq_along(segments)) {
        seqSetFilter(seqData, variant.sel=segments[i], verbose=FALSE)
        if (sum(seqGetFilter(seqData)$variant.sel) == 0) next
        iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
        assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
        files[i] <- tempfile()
        save(assoc, file=files[i])
    }
    files <- setdiff(files, "")

    assoc <- combineAssoc(files, "single")

    seqSetFilterChrom(seqData, include=1, verbose=FALSE)
    iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
    a <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    expect_equivalent(a, assoc)

    seqClose(seqData)
    unlink(files)
})


test_that("aggregate - GRangesList", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    data(segments)
    segfile <- .testSegFile(segments)
    segments <- segments[seqnames(segments) == 1]
    
    id <- which(seqGetData(seqData, "chromosome") == 1)
    pos <- seqGetData(seqData, "position")[id]
    bins <- cut(id, breaks=10)
    varList <- GRangesList(lapply(levels(bins), function(x) {
        ind <- which(bins == x)
        GRanges(seqnames=1, ranges=IRanges(start=pos[ind], width=1))
    }))
    files <- character(length(segments))
    for (i in seq_along(segments)) {
        vl <- subsetBySegment(varList, i, segfile)
        if (length(vl) == 0) next
        iterator <- SeqVarListIterator(seqData, vl, verbose=FALSE)
        assoc <- assocTestAggregate(iterator, nullmod, verbose=FALSE)
        files[i] <- tempfile()
        save(assoc, file=files[i])
        seqResetFilter(seqData, verbose=FALSE)
    }
    files <- setdiff(files, "")

    assoc <- combineAssoc(files, "aggregate")

    seqResetFilter(seqData, verbose=FALSE)
    iterator <- SeqVarListIterator(seqData, varList, verbose=FALSE)
    a <- assocTestAggregate(iterator, nullmod, verbose=FALSE)
    expect_equal(a, assoc)

    seqClose(seqData)
    unlink(files)
    unlink(segfile)
})

test_that("aggregate - GRanges", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    data(segments)
    segfile <- .testSegFile(segments)
    segments <- segments[seqnames(segments) == 1]
    
    id <- which(seqGetData(seqData, "chromosome") == 1)
    pos <- seqGetData(seqData, "position")[id]
    bins <- cut(id, breaks=10)
    varList <- do.call(c, (lapply(levels(bins), function(x) {
        ind <- which(bins == x)
        GRanges(seqnames=1, ranges=IRanges(start=pos[ind], width=1))
    })))
    files <- character(length(segments))
    for (i in seq_along(segments)) {
        vl <- subsetBySegment(varList, i, segfile)
        if (length(vl) == 0) next
        iterator <- SeqVarRangeIterator(seqData, vl, verbose=FALSE)
        assoc <- assocTestAggregate(iterator, nullmod, verbose=FALSE)
        files[i] <- tempfile()
        save(assoc, file=files[i])
        seqResetFilter(seqData, verbose=FALSE)
    }
    files <- setdiff(files, "")

    assoc <- combineAssoc(files, "aggregate")

    seqResetFilter(seqData, verbose=FALSE)
    iterator <- SeqVarRangeIterator(seqData, varList, verbose=FALSE)
    a <- assocTestAggregate(iterator, nullmod, verbose=FALSE)
    expect_equal(a, assoc)

    seqClose(seqData)
    unlink(files)
    unlink(segfile)
})


test_that("window", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)

    pos <- seqGetData(seqData, "position")[seqGetData(seqData, "chromosome") == 22]
    gr <- GRanges(seqnames=22, IRanges(start=pos, end=pos))
    size <- 1000
    shift <- 500
    segments <- GRanges(seqnames=22, ranges=IRanges(start=seq(1, max(pos), shift*1000), width=size*1.5*1000))
    segments <- subsetByOverlaps(segments, gr)
    
    segfile <- .testSegFile(segments)

    files <- character(length(segments))
    for (i in seq_along(segments)) {
        filterBySegment(seqData, i, segfile, pad.right=size*1000, verbose=FALSE)
        if (sum(seqGetFilter(seqData)$variant.sel) == 0) next
        iterator <- SeqVarWindowIterator(seqData, windowSize=size, windowShift=shift, verbose=FALSE)
        assoc <- assocTestAggregate(iterator, nullmod, verbose=FALSE)
        files[i] <- tempfile()
        save(assoc, file=files[i])
    }
    files <- setdiff(files, "")

    assoc <- combineAssoc(files, "window")

    seqSetFilterChrom(seqData, include=22, verbose=FALSE)
    iterator <- SeqVarWindowIterator(seqData, windowSize=size, windowShift=shift, verbose=FALSE)
    a <- assocTestAggregate(iterator, nullmod, verbose=FALSE)
    expect_equal(a, assoc)

    seqClose(seqData)
    unlink(files)
    unlink(segfile)
})


test_that("build 38", {
    data(chromosomes_hg38)
    n <- sample(1:23, 1)
    seg <- defineSegments(n=n, build="hg38")
    expect_equivalent(seg, chromosomes_hg38)
})
