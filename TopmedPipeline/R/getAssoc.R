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


formatAssocSingle <- function(seqData, assoc) {

    names(assoc)[names(assoc) %in% c("snpID", "variant.id")] <- "variantID"
    names(assoc) <- sub(".Stat", ".stat", names(assoc), fixed=TRUE)
    names(assoc) <- sub(".Pval", ".pval", names(assoc), fixed=TRUE)
    
    seqSetFilter(seqData, variant.id=assoc$variantID, verbose=FALSE)
    assoc$pos <- seqGetData(seqData, "position")
    if (!("chr" %in% names(assoc))) {
        assoc$chr <- seqGetData(seqData, "chromosome")
    }
    if (!("MAF" %in% names(assoc))) {
        assoc$MAF <- pmin(assoc$freq, 1 - assoc$freq)
        assoc$minor.allele <- ifelse(assoc$freq > 0.5, "ref", "alt")
    }
    
    init.cols <- c("variantID", "chr", "pos", "n", "MAF", "minor.allele")
    cols <- setdiff(names(assoc), c(init.cols, "freq"))
    assoc <- assoc[,c(init.cols, cols)]

    assoc
}
