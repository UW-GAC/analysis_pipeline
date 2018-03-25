#' Write a bed file to use as a track in a LocusZoom plot
#'
#' @param x data.frame with columns "chr", "start", "end"
#' @param file output file name
#' @param track.label character string to label BED track in plot
#'
#' @importFrom dplyr "%>%" mutate_ rename_ select_
#' @export
writeBED <- function(x, file, track.label="") {
    x <- x %>%
        mutate_(chr=~paste0("chr", chr)) %>%
        select_("chr", "start", "end")
    x$name <- track.label
    x$score <- ""
    x$strand <- ""
    x$thickStart <- ""
    x$thickEnd <- ""
    x$itemRgb <- "255,0,255"
    write.table(x, file=file, quote=FALSE, row.names=FALSE, col.names=FALSE ,sep="\t")
}


#' Write a metal-format file with variants to use in a LocusZoom plot
#'
#' @param x data.frame with columns "chr", "pos", "pval"
#' @param file output file name
#' 
#' @importFrom dplyr "%>%" mutate_ rename_ select_
#' @export
writeMETAL <- function(x, file) {
    x <- x %>%
        mutate_(MarkerName=~paste0("chr", chr, ":", pos)) %>%
        select_("MarkerName", "pval")
    names(x)[2] <- "P-value"

    write.table(x, file=file, quote=FALSE, row.names=FALSE, sep="\t")
}



#' Calculate LD
#'
#' @param gdsfile character string with GDS file name
#' @param variant.id vector of variant.id
#' @param ref.var variant.id of reference variant (if NULL, LD is calculated pairwise for all variants)
#' @param sample.id vector of sample.id to use (if NULL, all samples)
#' @return r^2
#'
#' @import SeqArray
#' @importFrom SeqVarTools altDosage
#' @importFrom stats cor
#' @export
calculateLD <- function(gdsfile, variant.id, ref.var=NULL, sample.id=NULL) {
    gds <- seqOpen(gdsfile)
    seqSetFilter(gds, variant.id=variant.id, sample.id=sample.id, verbose=FALSE)
    geno <- altDosage(gds)
    seqClose(gds)

    # correlations for just ref.var, or the entire matrix?
    if (is.null(ref.var)){
        geno.ref <- geno
    } else {
        geno.ref <- geno[,as.character(ref.var)]
    }

    # drop unnecessary dimensions (if ref.var is a single variant)
    r <- drop(cor(geno, geno.ref, use="pairwise.complete.obs"))
    r^2
}


#' Write a LD file to use in a LocusZoom plot
#'
#' @param x data.frame with columns "variant.id", "chr", "pos"
#' @param ld vector of LD
#' @param ref.var variant.id of the reference variant
#' @param file output file name
#' 
#' @importFrom dplyr "%>%" mutate_ select_
#' @export
writeLD <- function(x, ld, ref.var, file) {
    x <- x %>%
        mutate_(snp1=~paste0("chr", chr, ":", pos))
    x$snp2 <- x$snp1[x$variant.id == ref.var]
    x$dprime <- 0
    x$rsquare <- ld
    x <- x %>%
        select_("snp1", "snp2", "dprime", "rsquare")

    write.table(x, file=file, quote=FALSE, row.names=FALSE)
}
