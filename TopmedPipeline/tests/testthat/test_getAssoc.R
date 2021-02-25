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

.testAssoc <- function(n = 100, nchr = 23) {
  data.frame(
    chr = sort(rep_len(1:nchr, n)),
    pos = sample(1e8, n),
    start = sample(1e8, n),
    end = sample(1e8, n),
    stat = rnorm(n),
    pval = runif(n),
    MAC = sample(1:1000, n),
    freq = runif(n),
    stringsAsFactors = FALSE
  ) %>%
  mutate(chr = ifelse(chr == 23, "X", chr)) %>%
  mutate(chr = ordered(chr, levels = c(1:22, "X"))) %>%
  arrange(chr, pos) %>%
  mutate(id = 1:n()) %>%
  select(id, everything())
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
    expect_true("MAF" %in% names(assoc))
    idx <- a$freq < 0.5
    expect_equal(assoc$MAF[idx], a$freq[idx])
    expect_equal(assoc$MAF[!idx], 1 - a$freq[!idx])

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

  seqClose(seqData)
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
    seqClose(seqData)
    unlink(cvfile)
})


test_that("assocFilterByFile with one varfile", {
  filename <- tempfile()
  assoc <- .testAssoc()

  # All ids selected
  tmp <- assoc$id
  save(tmp, file=filename)
  out <- assocFilterByFile(assoc, filename)
  expect_equal(out$id, tmp)
  expect_equal(out, assoc)

  # Only keep ids on chromosome 1
  tmp <- assoc %>%
    filter(chr == 1) %>%
    pull(id)
    #tmp <- assoc$id
  save(tmp, file=filename)
  out <- assocFilterByFile(assoc, filename)
  expect_equal(out$id, tmp)

  # Assoc has fewer ids than in file.
  tmp <- assoc$id
  save(tmp, file = filename)
  assoc_subset <- assoc %>%
    filter(row_number() %in% 1:3)
  out <- assocFilterByFile(assoc_subset, filename)
  expect_equal(out$id, assoc_subset$id)

  # Assoc has repeated ids
  assoc_rep <- assoc %>%
    group_by(chr) %>%
    mutate(id = 1:n()) %>%
    ungroup()
  tmp <- assoc_rep$id[1:3]
  save(tmp, file = filename)
  chk <- assoc_rep %>%
    filter(id %in% tmp)
  out <- assocFilterByFile(assoc_rep, filename)
  expect_equal(out$id, chk$id)
  expect_equal(nrow(out), nrow(chk))
  expect_equal(out, chk)

  # No records in varfile
  tmp <- c()
  save(tmp, file = filename)
  out <- assocFilterByFile(assoc, filename)
  expect_equal(nrow(out), 0)

  # Different ids in varfile
  tmp <- max(assoc$id) + 1
  save(tmp, file = filename)
  out <- assocFilterByFile(assoc, filename)
  expect_equal(nrow(out), 0)

  unlink(filename)

})


test_that("assocFilterByFile with multiple varfiles", {
  file_pattern <- tempfile("assoc_chr _")
  assoc <- .testAssoc(nchr = 2, n = 6)

  tmp <- assoc %>%
    group_by(chr) %>%
    group_split()
  lapply(seq_along(tmp), function(x) {
    ids = tmp[[x]]$id
    chr = unique(tmp[[x]]$chr) %>% as.character()
    save(ids, file = insertChromString(file_pattern, chr))
  })
  out <- assocFilterByFile(assoc, file_pattern)
  expect_equal(out, assoc)

  # Extra var file
  extra_ids <- max(assoc$id) + 1
  extra_chr <- as.integer(as.character(max(assoc$chr))) + 1
  save(extra_ids, file=insertChromString(file_pattern, extra_chr))
  out <- assocFilterByFile(assoc, file_pattern)
  expect_equal(out$id, assoc$id)
  expect_equal(names(out), names(assoc))
  unlink(insertChromString(file_pattern, extra_chr))

  # Some ids selected.
  tmp <- assoc %>%
    group_by(chr) %>%
    dplyr::slice(1:2) %>%
    group_split()
  lapply(seq_along(tmp), function(x) {
    ids = tmp[[x]]$id
    chr = unique(tmp[[x]]$chr) %>% as.character()
    save(ids, file = insertChromString(file_pattern, chr))
  })

  out <- assocFilterByFile(assoc, file_pattern)
  expect_equal(out$id, tmp %>% bind_rows() %>% pull(id))
  expect_equal(names(out), names(assoc))

  # Only one chromosome in assoc but varfiles for other chromosomes exist.
  out <- assocFilterByFile(tmp[[1]], file_pattern)
  expect_equal(out$id, tmp[[1]]$id)
  expect_equal(names(out), names(assoc))

  # Assoc has fewer ids than in file across all chromosomes
  assoc_subset <- assoc %>%
    dplyr::group_by(chr) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  out <- assocFilterByFile(assoc_subset, file_pattern)
  expect_equal(out$id, assoc_subset$id)
  expect_equal(names(out), names(assoc))

  # Only chr1 varfile exists.
  unlink(insertChromString(file_pattern, 2))
  expect_warning(out <- assocFilterByFile(assoc, file_pattern), "missing varfile")
  expect_equal(out$id, assoc %>% filter(chr == 1) %>% dplyr::slice(1:2) %>% pull(id))
  expect_equal(names(out), names(assoc))
  #fail("Should there be a warning if a varfile is missing?")

  # No records in varfiles
  lapply(1:2, function(x) {
    ids = c()
    save(ids, file = insertChromString(file_pattern, x))
  })
  out <- assocFilterByFile(assoc, file_pattern)
  expect_equal(names(out), names(assoc))
  expect_equal(nrow(out), 0)

  # Chr1 file has chr2 ids and vice versa
  r = 1:2
  lapply(r, function(x) {
    ids <- tmp[[x]]$id
    save(ids, file = insertChromString(file_pattern, rev(r)[x]))
  })
  out <- assocFilterByFile(assoc, file_pattern)
  expect_equal(nrow(out), 0)

  # Different ids in varfiles
  tmp <- max(assoc$id) + c(1,2)
  lapply(seq_along(tmp), function(x) {
    ids <- tmp[x]
    save(ids, file = insertChromString(file_pattern, x))
  })
  out <- assocFilterByFile(assoc, file_pattern)
  expect_equal(names(out), names(assoc))
  expect_equal(nrow(out), 0)

  # Assoc has repeated ids
  assoc <- assoc %>%
    group_by(chr) %>%
    mutate(id = 1:n()) %>%
    ungroup()
  tmp <- assoc %>%
    group_by(chr) %>%
    dplyr::slice(1:2) %>%
    group_split()
  lapply(seq_along(tmp), function(x) {
    ids = tmp[[x]]$id
    chr = unique(tmp[[x]]$chr) %>% as.character()
    save(ids, file = insertChromString(file_pattern, chr))
  })

  out <- assocFilterByFile(assoc, file_pattern)
  expect_equal(names(out), names(assoc))
  expect_equal(out, tmp %>% bind_rows())

})
