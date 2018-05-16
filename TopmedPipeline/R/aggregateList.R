
#' @importFrom SeqVarTools expandedVariantIndex nAlleles variantInfo
#' @importFrom dplyr "%>%" rename_
.expandAlleles <- function(gds) {
    x <- variantInfo(gds, alleles=TRUE, expanded=TRUE) %>%
        rename_(allele="alt")
    x$nAlleles <- nAlleles(gds)[expandedVariantIndex(gds)]
    x
}


#' @importFrom dplyr "%>%" filter_ one_of select_
#' @importFrom stats setNames
.groupVariants <- function(variants, indexOnly) {

    ## columns to return
    extraCols <- if (indexOnly) character(0) else c("chr", "pos", "ref", "nAlleles", "allele")
    
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
#' @param variants A data.frame of variants with columns "group_id", "chr", "pos", "ref", "alt".
#' @param indexOnly Logical for whether to return only the "variant.id" and "allele.index" columns in the output (see Value).
#' @return A list of data frames, one for each group. A single variant may have multiple rows if multiple alternate alleles are selected. Each data frame contains the following columns:
#' \item{variant.id}{Unique identifier for the variant}
#' \item{chr}{Chromosome}
#' \item{pos}{Position in base pairs}
#' \item{ref}{Reference allele}
#' \item{nAlleles}{Total number of alleles for this variant}
#' \item{allele}{Alternate allele}
#' \item{allele.index}{Integer index of this allele (1=first alternate, 2=second alternate, etc.)}
#' @examples
#' library(SeqVarTools)
#' gds <- seqOpen(seqExampleFileName())
#' seqSetFilter(gds, variant.sel=seqGetData(gds, "chromosome") == 22)
#' variants <- data.frame(chr=seqGetData(gds, "chromosome"),
#'                        pos=seqGetData(gds, "position"),
#'                        ref=refChar(gds),
#'                        alt=altChar(gds, n=1),
#'                        stringsAsFactors=FALSE)
#' variants$group_id <- sample(LETTERS[1:2], nrow(variants), replace=TRUE)
#' aggregateListByAllele(gds, variants)
#' 
#' groups <- data.frame(group_id=LETTERS[1:2],
#'                      chr=22,
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
    stopifnot(all(c("group_id", "chr", "pos", "ref", "alt") %in% names(variants)))

    ## set filter to listed variants only
    filtOrig <- seqGetFilter(gds)
    gr <- GRanges(seqnames=variants$chr,
                  ranges=IRanges(variants$pos, variants$pos))
    seqSetFilter(gds, gr, verbose=FALSE)
  
    variants <- .expandAlleles(gds) %>%
        inner_join(variants, by=c("chr", "pos", "ref", allele="alt"))

    seqSetFilter(gds, sample.sel=filtOrig$sample.sel,
                 variant.sel=filtOrig$variant.sel, verbose=FALSE)

    .groupVariants(variants, indexOnly)
}


#' @param groups A data.frame of groups with column "group_id", "chr", "start", "end".
#' @rdname aggregateList
#'
#' @import SeqArray
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom dplyr "%>%" distinct_ left_join
#' @export
aggregateListByPosition <- function(gds, groups, indexOnly=FALSE) {
    stopifnot(all(c("group_id", "chr", "start", "end") %in% names(groups)))

    ## select only variants in requested regions
    filtOrig <- seqGetFilter(gds)
    gr <- GRanges(seqnames=groups$chr,
                  ranges=IRanges(groups$start, groups$end, names=groups$group_id))
    seqSetFilter(gds, gr, verbose=FALSE)

    variants <- .expandAlleles(gds)

    seqSetFilter(gds, sample.sel=filtOrig$sample.sel,
                 variant.sel=filtOrig$variant.sel, verbose=FALSE)
    
    ## find group_id for each variant
    vr <- GRanges(seqnames=variants$chr,
                  ranges=IRanges(variants$pos, variants$pos, names=variants$variant.id))
    ol <- findOverlaps(vr, gr)
    map <- data.frame(group_id=names(gr)[subjectHits(ol)],
                      variant.id=as.integer(names(vr))[queryHits(ol)],
                      stringsAsFactors=FALSE)

    variants <- distinct_(map) %>%
        left_join(variants, by="variant.id")
    .groupVariants(variants, indexOnly)
}



#' Aggregate variant lists
#'
#' Generate GRanges or GRangesList of variants for input to association tests
#'
#' These functions produce output suitable for defining a \code{\link{SeqVarRangeIterator}} or \code{\link{SeqVarListIterator}} object.
#'
#' @param variants A data.frame of variants with columns "group_id", "chr", "pos", "ref", "alt".
#' @return A GRangesList with one element per group
#' @examples
#' library(SeqVarTools)
#' gds <- seqOpen(seqExampleFileName())
#' seqSetFilter(gds, variant.sel=seqGetData(gds, "chromosome") == 22)
#' variants <- data.frame(chr=seqGetData(gds, "chromosome"),
#'                        pos=seqGetData(gds, "position"),
#'                        ref=refChar(gds),
#'                        alt=altChar(gds, n=1),
#'                        stringsAsFactors=FALSE)
#' variants$group_id <- sample(LETTERS[1:2], nrow(variants), replace=TRUE)
#' gr <- aggregateGRangesList(variants)
#' iterator <- SeqVarListIterator(gds, variantRanges=gr)
#' 
#' groups <- data.frame(group_id=LETTERS[1:2],
#'                      chr=22,
#'                      start=c(16000000, 2900000), 
#'                      end=c(30000000, 49000000),
#' 		     stringsAsFactors=FALSE)
#' gr <- aggregateGRanges(groups)
#' seqResetFilter(gds)
#' iterator <- SeqVarRangeIterator(gds, variantRanges=gr)
#' 
#' seqClose(gds)
#' @name aggregateGRanges
#'
#' @importFrom GenomicRanges GRanges GRangesList mcols<-
#' @importFrom IRanges IRanges
#' @importFrom stats setNames
#' @export
aggregateGRangesList <- function(variants) {
    stopifnot(all(c("group_id", "chr", "pos") %in% names(variants)))
    groups <- unique(variants$group_id)
    cols <- setdiff(names(variants), c("group_id", "chr", "pos"))
    GRangesList(lapply(setNames(groups, groups), function(g) {
        x <- variants[variants$group_id == g,]
        gr <- GRanges(seqnames=x$chr, ranges=IRanges(start=x$pos, width=1))
        mcols(gr) <- x[,cols]
        gr
    }))
}


#' @param groups A data.frame of groups with column "group_id", "chr", "start", "end".
#' @return A GRanges with one range per group
#' 
#' @rdname aggregateGRanges
#'
#' @importFrom GenomicRanges GRanges mcols<-
#' @importFrom IRanges IRanges
#' @export
aggregateGRanges <- function(groups) {
    stopifnot(all(c("group_id", "chr", "start", "end") %in% names(groups)))
    cols <- setdiff(names(groups), c("group_id", "chr", "start", "end"))
    gr <- GRanges(seqnames=groups$chr,
                  ranges=IRanges(start=groups$start, end=groups$end))
    names(gr) <- groups$group_id
    mcols(gr) <- groups[,cols]
    gr
}
