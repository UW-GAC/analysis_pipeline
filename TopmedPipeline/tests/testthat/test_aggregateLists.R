context("aggregateList tests")
library(SeqVarTools)
library(dplyr)

.opengds <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    seqOpen(seqExampleFileName())
}

.testVariants <- function(gds) {
    # triallelic snps are on chroms 21-22
    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=seqGetData(gds, "chromosome") %in% 21:22, verbose=FALSE)
    dat <- data.frame(chr=seqGetData(gds, "chromosome"),
                      pos=seqGetData(gds, "position"),
                      ref=refChar(gds),
                      alt=altChar(gds),
                      stringsAsFactors=FALSE)
    seqResetFilter(gds, verbose=FALSE)

    dat %>%
        mutate(alt=ifelse(nchar(alt) == 1, alt, substr(alt, 3, nchar(alt))),
               group_id=paste0("group", chr, sample(letters[1:2], n(), replace=TRUE)))
}

.testGroups <- function(gds) {
    dat <- .testVariants(gds) %>%
        group_by(group_id) %>%
        summarise(chr=unique(chr),
                  start=min(pos),
                  end=max(pos))

    # add overlapping groups
    dat2 <- dat %>%
        mutate(group_id=paste0(group_id, "c"),
               start=start + 100000,
               end=end + 100000)

    rbind(dat, dat2) %>%
        arrange(as.integer(chr), start)
}


test_that("expandAlleles", {
    gds <- .opengds()
    var <- .expandAlleles(gds)
    nAlt <- seqNumAllele(gds) - 1
    expect_equal(nrow(var), sum(nAlt))
    expect_equal(sum(var$allele.index == 2), sum(nAlt == 2))
    seqClose(gds)
})

test_that("aggregateListByAllele", {
    gds <- .opengds()
    variants <- .testVariants(gds)
    
    aggList <- aggregateListByAllele(gds, variants)
    expect_equal(unique(variants$group_id), names(aggList))
    
    df <- do.call(rbind, aggList)
    expect_true(all(c("variant.id", "allele.index") %in% names(df)))
    
    multi <- grepl(",", df$alt, fixed=TRUE)
    expect_true(all(df$allele.index[!multi] == 1))
    expect_true(all(df$allele.index[multi] == 2))
    
    seqClose(gds)
})

test_that("aggregateListByAllele returns correct columns", {
    gds <- .opengds()
    variants <- .testVariants(gds)
    
    aggList <- aggregateListByAllele(gds, variants)
    expect_equal(ncol(aggList[[1]]), 7)
    
    aggList <- aggregateListByAllele(gds, variants, indexOnly=TRUE)
    expect_equal(ncol(aggList[[1]]), 2)
    
    seqClose(gds)
})

test_that("aggregateListByAllele can handle multiple groups per variant", {
    gds <- .opengds()
    variants <- .testVariants(gds)
    groups <- unique(variants$group_id)
    variants <- filter_(variants, ~(group_id == groups[1])) %>%
        mutate_(group_id=~(groups[2])) %>%
        rbind(variants)
    
    aggList <- aggregateListByAllele(gds, variants)
    for (group in groups[1:2]) {
        g1 <- filter_(variants, ~(group_id == group)) %>%
            mutate(id=paste(chr, pos, ref))
        g2 <- aggList[[group]] %>%
            mutate_(id=~(paste(chr, pos, ref)))
        expect_true(setequal(g1$id, g2$id))
    }

    seqClose(gds)
})

test_that("aggregateListByPosition", {
    gds <- .opengds()
    groups <- .testGroups(gds)
    aggList <- aggregateListByPosition(gds, groups)
    expect_true(setequal(groups$group_id, names(aggList)))
    
    variants <- variantInfo(gds)   
    var.exp <- lapply(1:nrow(groups), function(i) {
        filter_(variants,
               ~(chr == groups$chr[i]),
               ~(pos >= groups$start[i]),
               ~(pos <= groups$end[i]))
    })
    names(var.exp) <- groups$group_id
    for (i in names(aggList)) {
        expect_true(setequal(var.exp[[i]]$variant.id, aggList[[i]]$variant.id))
    }
    
    seqClose(gds)
})

test_that("aggregateListByPosition gets all alternate alleles", {
    gds <- .opengds()
    groups <- .testGroups(gds)
    aggList <- aggregateListByPosition(gds, groups)

    lapply(aggList, function(x) {
        tmp <- filter_(x, ~(nAlleles == 3)) %>%
            group_by(variant.id) %>%
            summarise_(n=~(n()))
        expect_true(all(tmp$n == 2))
    })
    
    seqClose(gds)
})

test_that("aggregateGRangesList", {
    gds <- .opengds()
    variants <- .testVariants(gds)
    
    aggList <- aggregateGRangesList(variants)
    expect_equal(unique(variants$group_id), names(aggList))
    
    gr <- unlist(aggList)
    expect_true(all(c("ref", "alt") %in% names(GenomicRanges::mcols(gr))))
    
    expect_equal(length(gr), nrow(variants))
    
    seqClose(gds)
})

test_that("aggregateGRangesList can handle multiple groups per variant", {
    gds <- .opengds()
    variants <- .testVariants(gds)
    groups <- unique(variants$group_id)
    variants <- filter_(variants, ~(group_id == groups[1])) %>%
        mutate_(group_id=~(groups[2])) %>%
        rbind(variants)
    
    aggList <- aggregateGRangesList(variants)
    for (group in groups[1:2]) {
        g1 <- filter_(variants, ~(group_id == group)) %>%
            mutate(id=paste(chr, pos, ref))
        g2 <- aggList[[group]] %>%
            as.data.frame() %>%
            mutate_(id=~(paste(seqnames, start, ref)))
        expect_true(setequal(g1$id, g2$id))
    }

    seqClose(gds)
})

test_that("aggregateGRanges", {
    gds <- .opengds()
    groups <- .testGroups(gds)
    aggList <- aggregateGRanges(groups)
    expect_true(setequal(groups$group_id, names(aggList)))
    seqClose(gds)
})
