
.variantDF <- function(gds) {
    data.frame(variant.id=seqGetData(gds, "variant.id"),
               chromosome=seqGetData(gds, "chromosome"),
               position=seqGetData(gds, "position"),
               ref=refChar(gds),
               alt=altChar(gds),
               nAlleles=seqNumAllele(gds),
               stringsAsFactors=FALSE)
}


.expandAlleles <- function(gds) {
    variants <- .variantDF(gds)
    
    alleleCount <- sort(unique(variants$nAlleles))
    alleles <- lapply(alleleCount, function(n) {
        if (n == 2) {
            filter_(variants, ~(nAlleles == n)) %>%
                mutate_(allele.index=1, allele="alt") %>%
                select_("-alt")
        } else {
            alt.cols <- paste0("alt", 1:(n-1))
            filter_(variants, ~(nAlleles == n)) %>%
                separate_("alt", into=alt.cols, sep=",") %>%
                gather_("allele.index", "allele", alt.cols) %>%
                mutate_(allele.index=~(sub("alt", "", allele.index)))
        }
    })
    do.call(rbind, alleles) %>%
        arrange_("variant.id")
}


.groupVariants <- function(variants, indexOnly) {

    ## columns to return
    extraCols <- if (indexOnly) character(0) else c("chromosome", "position", "ref", "nAlleles", "allele")
    
    groups <- unique(variants$group_id)
    lapply(setNames(groups, groups), function(g) {
        filter_(variants, ~(group_id == g)) %>%
            select_(~(one_of(c("variant.id", extraCols, "allele.index"))))
    })
}


aggregateListByAllele <- function(gds, variants, indexOnly=FALSE) {
    stopifnot(all(c("group_id", "chromosome", "position", "ref", "alt") %in% names(variants)))

    ## set filter to listed variants only
    filtOrig <- seqGetFilter(gds)
    gr <- GRanges(seqnames=variants$chromosome,
                  ranges=IRanges(variants$position, variants$position))
    seqSetFilter(gds, gr, verbose=FALSE)
  
    variants <- .expandAlleles(gds) %>%
        inner_join(variants, by=c("chromosome", "position", "ref", allele="alt"))

    seqSetFilter(gds, sample.sel=filtOrig$sample.sel,
                 variant.sel=filtOrig$variant.sel, verbose=FALSE)

    .groupVariants(variants, indexOnly)
}


aggregateListByPosition <- function(gds, groups, indexOnly=FALSE) {
    stopifnot(all(c("group_id", "chromosome", "start", "end") %in% names(groups)))

    ## select only variants in requested regions
    filtOrig <- seqGetFilter(gds)
    gr <- GRanges(seqnames=groups$chromosome,
                  ranges=IRanges(groups$start, groups$end, names=groups$group_id))
    seqSetFilter(gds, gr, verbose=FALSE)

    variants <- .expandAlleles(gds)

    seqSetFilter(gds, sample.sel=filtOrig$sample.sel,
                 variant.sel=filtOrig$variant.sel, verbose=FALSE)
    
    ## find group_id for each variant
    vr <- GRanges(seqnames=variants$chromosome,
                  ranges=IRanges(variants$position, variants$position, names=variants$variant.id))
    ol <- findOverlaps(vr, gr)
    map <- data.frame(group_id=names(gr)[subjectHits(ol)],
                      variant.id=as.integer(names(vr))[queryHits(ol)],
                      stringsAsFactors=FALSE)

    variants <- distinct_(map) %>%
        left_join(variants, by="variant.id")
    .groupVariants(variants, indexOnly)
}
