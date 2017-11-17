context("getAssoc tests")
library(genesis2)
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
    fitNullModel2(sampleData(seqData), outcome="outcome", covars="sex", covMatList=grm, verbose=FALSE)
}
        

test_that("single related", {
    seqData <- SeqVarBlockIterator(.testData(), verbose=FALSE)
    nullmod <- .testNullModel(seqData, MM=TRUE)
    seqSetFilterChrom(seqData, include=1, verbose=FALSE)
    a1 <- assocTestMM2(seqData, nullmod, verbose=FALSE)
    seqSetFilterChrom(seqData, include=2, verbose=FALSE)
    a2 <- assocTestMM2(seqData, nullmod, verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])
    a <- rbind(a1, a2) %>%
        filter_(~!is.na(Wald.pval))

    assoc <- getAssoc(files, "single")
    expect_equal(as.character(assoc$chromosome), a$chromosome)
    expect_equal(assoc$position, a$position)
    expect_equal(assoc$start, assoc$pos)
    expect_equal(assoc$end, assoc$pos)
    expect_equal(assoc$stat, a$Wald.Stat)
    expect_equal(assoc$pval, a$Wald.pval)

    seqClose(seqData)
    unlink(files)
})


test_that("window, burden", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    seqSetFilterChrom(seqData, include=1, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a1 <- assocTestSeq2(seqData, nullmod, verbose=FALSE)
    seqSetFilterChrom(seqData, include=2, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a2 <- assocTestSeq2(seqData, nullmod, verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])
    a <- rbind(a1$results, a2$results) %>%
        filter_(~(n.site > 0))

    assoc <- getAssoc(files, "window")
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
    seqData <- SeqVarListIterator(seqData, variantRanges=agg, verbose=FALSE)
    a1 <- assocTestSeq2(seqData, nullmod, test="SKAT", verbose=FALSE)
    restoreFilter(seqData)
    agg <- GRangesList(gr[301:400], gr[401:500], gr[501:600])
    seqData <- SeqVarListIterator(seqData, variantRanges=agg, verbose=FALSE)
    a2 <- assocTestSeq2(seqData, nullmod, test="SKAT", verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])
    a <- rbind(a1$results, a2$results) %>%
        filter_(~(n.site > 0))

    assoc <- getAssoc(files, "aggregate")
    expect_true(setequal(assoc$pval, a$pval_0))
    expect_equal(as.character(assoc$chromosome[1]), a1$variantInfo[[1]]$chromosome[1])
    expect_true(assoc$position[1] > a1$variantInfo[[1]]$position[1] & assoc$position[1] < max(a1$variantInfo[[1]]$position))
    expect_equal(assoc$start[1], a1$variantInfo[[1]]$position[1])
    expect_equal(assoc$end[1], max(a1$variantInfo[[1]]$position))

    seqClose(seqData)
    unlink(files)
})


test_that("combine single", {
    seqData <- SeqVarBlockIterator(.testData(), verbose=FALSE)
    nullmod <- .testNullModel(seqData, MM=TRUE)
    seqSetFilterChrom(seqData, include=1, verbose=FALSE)
    a1 <- assocTestMM2(seqData, nullmod, verbose=FALSE)
    seqSetFilterChrom(seqData, include=2, verbose=FALSE)
    a2 <- assocTestMM2(seqData, nullmod, verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])

    assoc <- combineAssoc(files, "single")
    
    seqSetFilterChrom(seqData, include=1:2, verbose=FALSE)
    a <- assocTestMM2(seqData, nullmod, verbose=FALSE)
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
    a1 <- assocTestSeq2(seqData, nullmod, verbose=FALSE)
    seqSetFilterChrom(seqData, include=2, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a2 <- assocTestSeq2(seqData, nullmod, verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])

    assoc <- combineAssoc(files, "window")
    
    seqData <- .testData()
    seqSetFilterChrom(seqData, include=1:2, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a <- assocTestSeq2(seqData, nullmod, verbose=FALSE)
    # windows with n.site=0 are at end of combine
    #expect_equal(a, assoc)

    seqClose(seqData)
    unlink(files)
})

test_that("combine aggregate", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    gr <- granges(seqData)
    agg1 <- GRangesList(gr[1:100], gr[101:200], gr[201:300])
    seqData <- SeqVarListIterator(seqData, variantRanges=agg1, verbose=FALSE)
    a1 <- assocTestSeq2(seqData, nullmod, test="SKAT", verbose=FALSE)
    restoreFilter(seqData)
    agg2 <- GRangesList(gr[301:400], gr[401:500], gr[501:600])
    seqData <- SeqVarListIterator(seqData, variantRanges=agg2, verbose=FALSE)
    a2 <- assocTestSeq2(seqData, nullmod, test="SKAT", verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])

    assoc <- combineAssoc(files, "aggregate")
    
    restoreFilter(seqData)
    seqData <- SeqVarListIterator(seqData, variantRanges=c(agg1, agg2), verbose=FALSE)
    a <- assocTestSeq2(seqData, nullmod, test="SKAT", verbose=FALSE)
    expect_equal(a, assoc)

    seqClose(seqData)
    unlink(files)
})

test_that("arrange_chr_pos", {
    x <- data.frame(chromosome=c(10,2,"X",1,1), position=c(1,1,1,2,1), stringsAsFactors=FALSE)
    xa <- x[c(5,4,2,1,3),]; rownames(xa) <- 1:5
    expect_equal(.arrange_chr_pos(x), xa)
    names(x) <- names(xa) <- c("a", "b")
    expect_equal(.arrange_chr_pos(x, chr="a", pos="b"), xa)
})

test_that("index_chr_pos", {
    x <- list(data.frame(chromosome=2, position=1:5),
              data.frame(chromosome=1, position=6:10),
              data.frame(chromosome="X", position=1:5, stringsAsFactors=FALSE),
              data.frame(chromosome=1, position=1:5))
    expect_equal(.index_chr_pos(x), c(4,2,1,3))
    xa <- lapply(x, function(xx) {names(xx) <- c("a", "b"); xx})
    expect_equal(.index_chr_pos(xa, chr="a", pos="b"), c(4,2,1,3))
})
    
test_that("combine out of order", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData, MM=TRUE)

    seqSetFilter(seqData, variant.sel=1:10, verbose=FALSE)
    seqData <- SeqVarBlockIterator(seqData, verbose=FALSE)
    a1 <- assocTestMM2(seqData, nullmod, verbose=FALSE)
    seqSetFilter(seqData, variant.sel=11:20, verbose=FALSE)
    seqData <- SeqVarBlockIterator(seqData, verbose=FALSE)
    a2 <- assocTestMM2(seqData, nullmod, verbose=FALSE)
    files <- file.path(tempdir(), c("a", "b"))
    save(a1, file=files[1])
    save(a2, file=files[2])
    assoc <- combineAssoc(files, "single")
    a <- combineAssoc(files[c(2,1)], "single")
    expect_equal(a, assoc)

    seqSetFilterChrom(seqData, include=1, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a1 <- assocTestSeq2(seqData, nullmod, verbose=FALSE)
    seqSetFilterChrom(seqData, include=2, verbose=FALSE)
    seqData <- SeqVarWindowIterator(seqData, verbose=FALSE)
    a2 <- assocTestSeq2(seqData, nullmod, verbose=FALSE)
    save(a1, file=files[1])
    save(a2, file=files[2])
    assoc <- combineAssoc(files, "window")
    a <- combineAssoc(files[c(2,1)], "window") 
    #expect_equal(a, assoc)

    seqResetFilter(seqData, verbose=FALSE)
    gr <- granges(seqData)
    agg1 <- GRangesList(gr[1:100], gr[101:200], gr[201:300])
    seqData <- SeqVarListIterator(seqData, variantRanges=agg1, verbose=FALSE)
    a1 <- assocTestSeq2(seqData, nullmod, test="SKAT", verbose=FALSE)
    restoreFilter(seqData)
    agg2 <- GRangesList(gr[301:400], gr[401:500], gr[501:600])
    seqData <- SeqVarListIterator(seqData, variantRanges=agg2, verbose=FALSE)
    a2 <- assocTestSeq2(seqData, nullmod, test="SKAT", verbose=FALSE)
    save(a1, file=files[1])
    save(a2, file=files[2])
    assoc <- combineAssoc(files, "aggregate")
    agg <- c(agg1, agg2)[c(2,1,6,4,5,3)]
    restoreFilter(seqData)
    seqData <- SeqVarListIterator(seqData, variantRanges=agg, verbose=FALSE)
    a <- assocTestSeq2(seqData, nullmod, test="SKAT", verbose=FALSE)
    #expect_equal(a, assoc)
    
    seqClose(seqData)
    unlink(files)
})


test_that("omitKnownHits", {
    assoc <- data.frame(chromosome=c(rep(1,100), rep(2,100)),
                        position=c(sort(sample(1:10000, 100)), sort(sample(1:10000, 100))))
    hits <- assoc[sample(1:nrow(assoc), 20),]
    res <- omitKnownHits(assoc, hits, flank=0)
    .pos <- function(x) paste0(x$chr, x$pos)
    expect_equal(setdiff(.pos(assoc), .pos(hits)), .pos(res))
    
    res <- omitKnownHits(assoc, hits, flank=1)
    expect_true(nrow(res) < nrow(assoc))
})
