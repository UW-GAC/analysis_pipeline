context("filterVariants tests")
library(dplyr)
library(gdsfmt)
library(GenomicRanges)

.testData <- function() {
    showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- seqExampleFileName("gds")
    seqOpen(gdsfile)
}

.testBinaryData <- function() {
    gds <- .testData()
    samp <- seqGetData(gds, "sample.id")
    annot <- AnnotatedDataFrame(data.frame(sample.id=samp,
                                           outcome=rbinom(length(samp), 1, 0.3),
                                           stringsAsFactors=FALSE))
    SeqVarData(gds, annot)
}

.testSegFile <- function() {
    data(segments)
    seg.df <- as.data.frame(segments) %>%
        dplyr::rename(chromosome=seqnames) %>%
        dplyr::select(chromosome, start, end)
    segfile <- tempfile()
    write.table(seg.df, file=segfile, quote=FALSE, sep="\t", row.names=FALSE)
    segfile
}

.testVarList <- function(gds) {
    var.id <- seqGetData(gds, "variant.id")
    chr <- seqGetData(gds, "chromosome")
    pos <- seqGetData(gds, "position")
    ind1 <- which(chr == 1)[1:10]
    ind2 <- which(chr == 1)[11:20]
    ind3 <- which(chr == 2)[1:10]
    lapply(list(ind1, ind2, ind3), function(x) {
        data.frame(variant.id=var.id[x], chromosome=chr[x], position=pos[x], allele.index=1, stringsAsFactors=FALSE)
    })
}

test_that("getSegments", {
    segfile <- .testSegFile()
    seg2 <- getSegments(segfile)
    expect_equivalent(seg2, segments)
    unlink(segfile)
})

test_that("subsetBySegment", {
    gds <- .testData()
    varList <- .testVarList(gds)
    segfile <- .testSegFile()
    segments <- getSegments(segfile)
    exp <- sapply(varList, function(x) {
        chr <- x$chromosome[1]
        pos <- x$position[1]
        seg.chr <- as.character(seqnames(segments[1]))
        seg.start <- as.integer(BiocGenerics::start(segments[1]))
        seg.end <- as.integer(BiocGenerics::end(segments[1]))
        chr == seg.chr & pos >= seg.start & pos <= seg.end
    })
    ss <- subsetBySegment(varList, 1, segfile)
    expect_equivalent(varList[exp], ss)
    
    seqClose(gds)
    unlink(segfile)
})

test_that("filterBySegment", {
    gds <- .testData()
    segfile <- .testSegFile()
    segments <- getSegments(segfile)
    gr <- granges(gds)
    ol <- findOverlaps(gr, segments[1])
    filterBySegment(gds, 1, segfile, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), length(queryHits(ol)))

    seqClose(gds)
    unlink(segfile)
})


test_that("filterByFile", {
    gds <- .testData()
    id <- seqGetData(gds, "variant.id")[1:100]
    idfile <- tempfile()
    save(id, file=idfile)
    filterByFile(gds, idfile, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), 100)

    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=1:10, verbose=FALSE)
    filterByFile(gds, idfile, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), 10)

    seqClose(gds)
    unlink(idfile)
})


test_that("filterByChrom", {
    gds <- .testData()
    chr <- seqGetData(gds, "chromosome")
    filterByChrom(gds, 1, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), sum(chr == 1))

    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=1:10, verbose=FALSE)
    filterByChrom(gds, 1, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), 10)

    seqClose(gds)
})


test_that("filterByPass", {
    gds <- .testData()
    filt <- seqGetData(gds, "annotation/filter")
    filterByPass(gds, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), sum(filt == "PASS"))

    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=1:10, verbose=FALSE)
    filterByPass(gds, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), 10)

    seqClose(gds)
})


test_that("filterBySNV", {
    gds <- .testData()
    snv <- isSNV(gds, biallelic=TRUE)
    filterBySNV(gds, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), sum(snv))

    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=1:10, verbose=FALSE)
    filterBySNV(gds, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), 10)

    seqClose(gds)
})


test_that("filterByMAF", {
    gds <- .testData()
    freq <- seqAlleleFreq(gds)
    maf <- pmin(freq, 1-freq)
    filterByMAF(gds, maf.min=0.1, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), sum(maf >= 0.1))

    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=1:10, verbose=FALSE)
    filterByMAF(gds, maf.min=0.1, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), sum(maf[1:10] >= 0.1))

    seqClose(gds)
})


test_that("filterByMAC", {
    gds <- .testData()
    cnt1 <- alleleCount(gds)
    cnt2 <- alleleCount(gds, n=1)
    cnt <- round(pmin(cnt1, cnt2))
    x <- .calcMAC(gds, sample.id=NULL)
    expect_equal(x, cnt)
    filterByMAC(gds, mac.min=5, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), sum(cnt >= 5))

    seqClose(gds)
})


test_that("minAltFreq", {
    f <- list(c(1,0.5,0.2), c(1,0,0), c(1,0.1,0,1))
    expect_equal(.minAltFreq(f), c(0.2,NA,0.1))
})

test_that("filterByRare", {
    gds <- .testData()
    freq <- seqAlleleFreq(gds, ref.allele=1)
    filterByRare(gds, af.max=0.1, verbose=FALSE)
    expect_equal(sum(seqGetFilter(gds)$variant.sel), sum(freq > 0 & freq <= 0.1))
    seqClose(gds)
})



test_that("filterByPCAcorr", {
    gds <- .testData()
    filt <- get(data(pcaSnpFilters.hg19, package="GWASTools"))
    filt.gr <- GRanges(seqnames=filt$chrom, ranges=IRanges(start=filt$start.base, end=filt$end.base))
    exp <- GenomicRanges::setdiff(granges(gds), filt.gr)

    filterByPCAcorr(gds, build="hg19", verbose=FALSE)
    expect_equal(length(GenomicRanges::setdiff(granges(gds), exp)), 0)
    
    seqClose(gds)
})
