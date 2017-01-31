context("getAssoc tests")
library(GENESIS)
library(gdsfmt)
library(SeqVarTools)
library(Biobase)

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
        fitNullMM(sampleData(seqData), outcome="outcome", covars="sex", covMatList=grm, verbose=FALSE)
    } else {
        fitNullReg(sampleData(seqData), outcome="outcome", covars="sex", verbose=FALSE)
    }
}

test_that("formatAssocSingle related", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData, MM=TRUE)
    assoc <- assocTestMM(seqData, nullmod, verbose=FALSE)
    assoc <- formatAssocSingle(seqData, assoc)
    seqClose(seqData)
    expect_true(all(c("variantID", "chr", "pos", "MAF") %in% names(assoc)))
})
   
test_that("formatAssocSingle unrelated", {
    seqData <- .testData()
    seqSetFilter(seqData, variant.sel=1:10, verbose=FALSE)
    assoc <- regression(seqData, outcome="outcome", covar="sex")
    assoc <- formatAssocSingle(seqData, assoc)
    seqClose(seqData)
    expect_true(all(c("variantID", "chr", "pos", "MAF") %in% names(assoc)))
})
                    
    


test_that("single unrelated", {
    seqData <- .testData()
    chr <- seqGetData(seqData, "chromosome")
    seqSetFilter(seqData, variant.sel=(chr == 1), verbose=FALSE)
    a1 <- formatAssocSingle(seqData, regression(seqData, outcome="outcome"))
    seqSetFilter(seqData, variant.sel=(chr == 2), verbose=FALSE)
    a2 <- formatAssocSingle(seqData, regression(seqData, outcome="outcome"))
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])
    a <- rbind(a1, a2) %>%
        filter_(~!is.na(Wald.pval))

    assoc <- getAssoc(files, "single")
    expect_equal(as.character(assoc$chr), a$chr)
    expect_equal(assoc$pos, a$pos)
    expect_equal(assoc$stat, a$Wald.stat)
    expect_equal(assoc$pval, a$Wald.pval)

    seqClose(seqData)
    unlink(files)
})


test_that("single related", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData, MM=TRUE)
    a1 <- formatAssocSingle(seqData, assocTestMM(seqData, nullmod, chromosome=1, verbose=FALSE))
    a2 <- formatAssocSingle(seqData, assocTestMM(seqData, nullmod, chromosome=2, verbose=FALSE))
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])
    a <- rbind(a1, a2) %>%
        filter_(~!is.na(Wald.pval))

    assoc <- getAssoc(files, "single")
    expect_equal(as.character(assoc$chr), a$chr)
    expect_equal(assoc$pos, a$pos)
    expect_equal(assoc$stat, a$Wald.stat)
    expect_equal(assoc$pval, a$Wald.pval)

    seqClose(seqData)
    unlink(files)
})


test_that("window, burden", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    a1 <- assocTestSeqWindow(seqData, nullmod, chromosome=1, verbose=FALSE)
    a2 <- assocTestSeqWindow(seqData, nullmod, chromosome=2, verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])
    a <- rbind(a1$results, a2$results) %>%
        filter_(~(n.site > 0), ~(dup == 0))

    assoc <- getAssoc(files, "window")
    expect_equal(as.character(assoc$chr), a$chr)
    expect_true(all(assoc$pos > a$window.start & assoc$pos < a$window.stop))
    expect_equal(assoc$stat, a$Score.stat)
    expect_equal(assoc$pval, a$Score.pval)

    seqClose(seqData)
    unlink(files)
})


test_that("aggregate, skat", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    agg <- list(data.frame(variant.id=1:100, allele.index=1),
                 data.frame(variant.id=101:200, allele.index=1),
                 data.frame(variant.id=201:300, allele.index=1))
    a1 <- assocTestSeq(seqData, nullmod, agg, test="SKAT", verbose=FALSE)
    agg <- list(data.frame(variant.id=301:400, allele.index=1),
                 data.frame(variant.id=401:500, allele.index=1),
                 data.frame(variant.id=501:600, allele.index=1))
    a2 <- assocTestSeq(seqData, nullmod, agg, test="SKAT", verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])
    a <- rbind(a1$results, a2$results) %>%
        filter_(~(n.site > 0))

    assoc <- getAssoc(files, "aggregate")
    expect_equal(assoc$pval, a$pval_0)
    expect_equal(as.character(assoc$chr[1]), a1$variantInfo[[1]]$chr[1])
    expect_equal(assoc$pos[1], a1$variantInfo[[1]]$pos[1])

    seqClose(seqData)
    unlink(files)
})


test_that("combine single", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData, MM=TRUE)
    a1 <- formatAssocSingle(seqData, assocTestMM(seqData, nullmod, chromosome=1, verbose=FALSE))
    a2 <- formatAssocSingle(seqData, assocTestMM(seqData, nullmod, chromosome=2, verbose=FALSE))
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])

    assoc <- combineAssoc(files, "single")
    
    a <- formatAssocSingle(seqData, assocTestMM(seqData, nullmod, chromosome=1:2, verbose=FALSE))
    expect_equal(a, assoc)

    seqClose(seqData)
    unlink(files)
})

test_that("combine window", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    a1 <- assocTestSeqWindow(seqData, nullmod, chromosome=1, verbose=FALSE)
    a2 <- assocTestSeqWindow(seqData, nullmod, chromosome=2, verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])

    assoc <- combineAssoc(files, "window")
    
    a <- assocTestSeqWindow(seqData, nullmod, chromosome=1:2, verbose=FALSE)
    expect_equal(a, assoc)

    seqClose(seqData)
    unlink(files)
})

test_that("combine aggregate", {
    seqData <- .testData()
    nullmod <- .testNullModel(seqData)
    agg1 <- list(data.frame(variant.id=1:100, allele.index=1),
                 data.frame(variant.id=101:200, allele.index=1),
                 data.frame(variant.id=201:300, allele.index=1))
    a1 <- assocTestSeq(seqData, nullmod, agg1, test="SKAT", verbose=FALSE)
    agg2 <- list(data.frame(variant.id=301:400, allele.index=1),
                 data.frame(variant.id=401:500, allele.index=1),
                 data.frame(variant.id=501:600, allele.index=1))
    a2 <- assocTestSeq(seqData, nullmod, agg2, test="SKAT", verbose=FALSE)
    files <- c(tempfile(), tempfile())
    save(a1, file=files[1])
    save(a2, file=files[2])

    assoc <- combineAssoc(files, "aggregate")

    a <- assocTestSeq(seqData, nullmod, c(agg1, agg2), test="SKAT", verbose=FALSE)
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
