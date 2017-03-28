#' Combine association test results
#'
#' Combine association test results from multiple files into a single object.
#' Useful for combining per-segment results into a single file per chromosome.
#' 
#' @param files Vector of file names with association test results
#' @param assoc_type Type of association test ("single", "aggregate", "window")
#' @return Association test object
#'
#' @importFrom dplyr "%>%" distinct_ filter_ group_by_
#' @export
combineAssoc <- function(files, assoc_type) {
    stopifnot(assoc_type %in% c("single", "aggregate", "window"))
    x <- lapply(unname(files), getobj)
    if (assoc_type == "single") {
        assoc <- do.call(rbind, x)
    } else if (assoc_type  == "aggregate") {
        assoc <- x[[1]][c("param", "nsample")]
        assoc$results <- do.call(rbind, lapply(x, function(y) y$results))
        assoc$variantInfo <- do.call(c, lapply(x, function(y) y$variantInfo))
    } else if (assoc_type == "window") {
        assoc <- x[[1]][c("param", "window", "nsample")]
        for (v in c("results", "variantInfo")) {
            assoc[[v]] <- do.call(rbind, lapply(x, function(y) y[[v]])) %>%
                distinct_()
        }
        assoc$results <- assoc$results %>%
            group_by_("chr", "window.start", "window.stop") %>%
            filter_(~(n.site == max(n.site)), ~(!duplicated(n.site))) %>%
            as.data.frame()
    }
    assoc
}

#' Get association test results
#'
#' Return association test results in a standard format
#'
#' Read association test results in multiple files and combine all into a single
#' data frame with standard column names.
#' 
#' @inheritParams combineAssoc
#' @return data.frame including standard columns ("chr", "pos", "stat", "pval")
#'
#' @importFrom dplyr "%>%" filter_ left_join mutate_ n rename_ select_
#' @export
getAssoc <- function(files, assoc_type) {
    stopifnot(assoc_type %in% c("single", "aggregate", "window"))
    assoc <- do.call(rbind, lapply(unname(files), function(f) {
        x <- getobj(f)
        if (assoc_type  == "aggregate") {
            tmp <- x$results %>%
                mutate_(group_id=~(1:n())) %>%
                filter_(~(n.site > 0))
            group.info <- do.call(rbind, lapply(tmp$group_id, function(g) {
                grp <- x$variantInfo[[g]][1, c("chr", "pos"), drop=FALSE]
                grp$group_id <- g
                grp
            }))
            x <- left_join(tmp, group.info, by="group_id")
        } else if (assoc_type == "window") {
            x <- filter_(x$results, ~(n.site > 0), ~(dup == 0)) %>%
                mutate_(pos=~(floor((window.start + window.stop)/2)))
        }
        x
    }))
    
    if ("pval_0" %in% names(assoc)) {
        ## SKAT
        pval.col <- if ("pval_SKATO" %in% names(assoc)) "pval_SKATO" else "pval_0"
        assoc <- select_(assoc, "chr", "pos", pval.col) %>%
            rename_(pval=pval.col)
    } else {
        ## burden or single
        assoc <- select_(assoc, "chr", "pos", ~ends_with("stat"), ~ends_with("pval"))
        names(assoc)[3:4] <- c("stat", "pval")
    }
    assoc <- filter_(assoc, ~(!is.na(pval))) %>%
        mutate_(chr=~factor(chr, levels=c(1:22, "X")))
    assoc
}

#' Format single-variant assocation test results
#'
#' Return association test results for single-variant tests in a standard format
#'
#' Ensures that single-variant association test results from different sources
#' (\code{\link[GENESIS]{assocTestMM}}, \code{\link[SeqVarTools]{regression}})
#' will all have a standard format.
#'
#' @param seqData A \code{\link[SeqArray]{SeqVarGDSClass}} object
#'   (needed to get chromosome and position if not present)
#' @param assoc data.frame with assocation test results
#' @return data.frame including standard columns ("variantID", "chr", "pos", "n", "MAF", "minor.allele")
#' 
#' @import SeqArray
#' @export
formatAssocSingle <- function(seqData, assoc) {

    names(assoc)[names(assoc) %in% c("snpID", "variant.id")] <- "variantID"
    names(assoc) <- sub(".Stat", ".stat", names(assoc), fixed=TRUE)
    names(assoc) <- sub(".Pval", ".pval", names(assoc), fixed=TRUE)
    
    seqSetFilter(seqData, variant.id=assoc$variantID, action="push+set", verbose=FALSE)
    assoc$pos <- seqGetData(seqData, "position")
    if (!("chr" %in% names(assoc))) {
        assoc$chr <- seqGetData(seqData, "chromosome")
    }
    if (!("MAF" %in% names(assoc))) {
        assoc$MAF <- pmin(assoc$freq, 1 - assoc$freq)
        assoc$minor.allele <- ifelse(assoc$freq > 0.5, "ref", "alt")
    }
    seqSetFilter(seqData, action="pop", verbose=FALSE)
    
    init.cols <- c("variantID", "chr", "pos", "n", "MAF", "minor.allele")
    cols <- setdiff(names(assoc), c(init.cols, "freq"))
    assoc <- assoc[,c(init.cols, cols)]

    assoc
}


#' Omit known hits from an association test data frame
#'
#' @param assoc data.frame with assocation test results (including columns chr, pos)
#' @param hits data.frame with known hits (including columns chr, pos)
#' @param flank Number of kb on either side of each known hit to exclude
#' @return data.frame with assocation test results not in regions around known hits
#'
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
#' @export
omitKnownHits <- function(assoc, hits, flank=500) {
    stopifnot(all(c("chr", "pos") %in% names(hits)))
    assoc.gr <- GRanges(seqnames=assoc$chr, ranges=IRanges(start=assoc$pos, end=assoc$pos))
    hits.gr <- GRanges(seqnames=hits$chr, ranges=IRanges(start=hits$pos-(flank*1000),
                                              end=hits$pos+(flank*1000)))
    ol <- findOverlaps(assoc.gr, hits.gr)
    assoc[-queryHits(ol),]
}
