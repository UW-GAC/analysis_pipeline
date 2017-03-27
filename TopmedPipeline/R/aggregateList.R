#' @import SeqArray
#' @importFrom SeqVarTools refChar altChar
.variantDF <- function(gds) {
    data.frame(variant.id=seqGetData(gds, "variant.id"),
               chromosome=seqGetData(gds, "chromosome"),
               position=seqGetData(gds, "position"),
               ref=refChar(gds),
               alt=altChar(gds),
               nAlleles=seqNumAllele(gds),
               stringsAsFactors=FALSE)
}


#' @importFrom dplyr "%>%" arrange_ filter_ mutate_ select_
#' @importFrom tidyr gather_ separate_ 
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


#' @importFrom dplyr "%>%" filter_ one_of select_
#' @importFrom stats setNames
.groupVariants <- function(variants, indexOnly) {

    ## columns to return
    extraCols <- if (indexOnly) character(0) else c("chromosome", "position", "ref", "nAlleles", "allele")
    
    groups <- unique(variants$group_id)
    lapply(setNames(groups, groups), function(g) {
        filter_(variants, ~(group_id == g)) %>%
            select_(~(one_of(c("variant.id", extraCols, "allele.index"))))
    })
}

#' Aggregate variant lists
#'
#' Generate lists of variants for input to association tests
#'
#' These functions produce output suitable for providing to \code{\link[GENESIS]{assocTestSeq}} in the \pkg{\link[GENESIS]{GENESIS}} package.
#'
#' @param gds A \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @param variants A data.frame of variants with columns "group_id", "chromosome", "position", "ref", "alt".
#' @param indexOnly Logical for whether to return only the "variant.id" and "allele.index" columns in the output (see Value).
#' @return A list of data frames, one for each group. A single variant may have multiple rows if multiple alternate alleles are selected. Each data frame contains the following columns:
#' \item{variant.id}{Unique identifier for the variant}
#' \item{chromosome}{Chromosome}
#' \item{position}{Position in base pairs}
#' \item{ref}{Reference allele}
#' \item{nAlleles}{Total number of alleles for this variant}
#' \item{allele}{Alternate allele}
#' \item{allele.index}{Integer index of this allele (1=first alternate, 2=second alternate, etc.)}
#' @examples
#' library(SeqVarTools)
#' gds <- seqOpen(seqExampleFileName())
#' seqSetFilter(gds, variant.sel=seqGetData(gds, "chromosome") == 22)
#' variants <- data.frame(chromosome=seqGetData(gds, "chromosome"),
#'                        position=seqGetData(gds, "position"),
#'                        ref=refChar(gds),
#'                        alt=altChar(gds, n=1),
#'                        stringsAsFactors=FALSE)
#' variants$group_id <- sample(LETTERS[1:2], nrow(variants), replace=TRUE)
#' aggregateListByAllele(gds, variants)
#' 
#' groups <- data.frame(group_id=LETTERS[1:2],
#'                      chromosome=22,
#'                      start=c(16000000, 2900000), 
#'                      end=c(30000000, 49000000),
#' 		     stringsAsFactors=FALSE)
#' aggregateListByPosition(gds, groups)
#' 
#' seqClose(gds)
#' @name aggregateList
#'
#' @import SeqArray
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom dplyr "%>%" inner_join
#' @export
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


#' @param groups A data.frame of groups with column "group_id", "chromosome", "start", "end".
#' @rdname aggregateList
#'
#' @import SeqArray
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom dplyr "%>%" distinct_ left_join
#' @export
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
