context("aggregateList tests")
library(dplyr)

.testVariants <- function(gds) {
    # triallelic snps are on chroms 21-22
    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=seqGetData(gds, "chromosome") %in% 21:22, verbose=FALSE)
    dat <- data.frame(chromosome=seqGetData(gds, "chromosome"),
                      position=seqGetData(gds, "position"),
                      ref=refChar(gds),
                      alt=altChar(gds),
                      stringsAsFactors=FALSE)
    seqResetFilter(gds, verbose=FALSE)

    dat %>%
        mutate(alt=ifelse(nchar(alt) == 1, alt, substr(alt, 3, nchar(alt))),
               group_id=paste0("group", chromosome, sample(letters[1:2], n(), replace=TRUE)))
}

.testGroups <- function(gds) {
    dat <- .testVariants(gds) %>%
        group_by(group_id) %>%
        summarise(chromosome=unique(chromosome),
                  start=min(position),
                  end=max(position))

    # add overlapping groups
    dat2 <- dat %>%
        mutate(group_id=paste0(group_id, "c"),
               start=start + 100000,
               end=end + 100000)

    rbind(dat, dat2) %>%
        arrange(as.integer(chromosome), start)
}


test_that("expandAlleles", {
    gds <- seqOpen(seqExampleFileName())
    var <- .expandAlleles(gds)
    nAlt <- seqNumAllele(gds) - 1
    expect_equal(nrow(var), sum(nAlt))
    expect_equal(sum(var$allele.index == 2), sum(nAlt == 2))
    seqClose(gds)
})

test_that("aggregateListByAllele", {
    gds <- seqOpen(seqExampleFileName())
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
    gds <- seqOpen(seqExampleFileName())
    variants <- .testVariants(gds)
    
    aggList <- aggregateListByAllele(gds, variants)
    expect_equal(ncol(aggList[[1]]), 7)
    
    aggList <- aggregateListByAllele(gds, variants, indexOnly=TRUE)
    expect_equal(ncol(aggList[[1]]), 2)
    
    seqClose(gds)
})

test_that("aggregateListByAllele can handle multiple groups per variant", {
    gds <- seqOpen(seqExampleFileName())
    variants <- .testVariants(gds)
    groups <- unique(variants$group_id)
    variants <- filter(variants, group_id == groups[1]) %>%
        mutate(group_id = groups[2]) %>%
        rbind(variants)
    
    aggList <- aggregateListByAllele(gds, variants)
    for (group in groups[1:2]) {
        g1 <- filter(variants, group_id == group) %>%
            mutate(id=paste(chromosome, position, ref))
        g2 <- aggList[[group]] %>%
            mutate(id=paste(chromosome, position, ref))
        expect_true(setequal(g1$id, g2$id))
    }

    seqClose(gds)
})

test_that("aggregateListByPosition", {
    gds <- seqOpen(seqExampleFileName())
    groups <- .testGroups(gds)
    aggList <- aggregateListByPosition(gds, groups)
    expect_true(setequal(groups$group_id, names(aggList)))
    
    variants <- .variantDF(gds)   
    var.exp <- lapply(1:nrow(groups), function(i) {
        filter(variants,
               chromosome == groups$chromosome[i],
               position >= groups$start[i],
               position <= groups$end[i])
    })
    names(var.exp) <- groups$group_id
    for (i in names(aggList)) {
        expect_true(setequal(var.exp[[i]]$variant.id, aggList[[i]]$variant.id))
    }
    
    seqClose(gds)
})

test_that("aggregateListByPosition gets all alternate alleles", {
    gds <- seqOpen(seqExampleFileName())
    groups <- .testGroups(gds)
    aggList <- aggregateListByPosition(gds, groups)

    lapply(aggList, function(x) {
        tmp <- filter(x, nAlleles == 3) %>%
            group_by(variant.id) %>%
            summarise(n=n())
        expect_true(all(tmp$n == 2))
    })
    
    seqClose(gds)
})
