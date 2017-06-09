#' Write a bed file to use as a track in a LocusZoom plot
#'
#' @param x data.frame with columns "chr", "start", "end"
#' @param file output file name
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
