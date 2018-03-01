#' Add variantData to a SeqVarData object
#'
#' @param gds A \code{\link{SeqVarData}} object
#' @param variants A data.frame with columns "variant.id", or ("chr", "pos", "ref", "alt")
#' @return \code{gds} with data from \code{variants} in the \code{variantData} slot
#'
#' @import SeqVarTools
#' @importFrom dplyr "%>%" left_join
#' @export
addVariantData <- function(gds, variants) {
    if (!("variant.id" %in% names(variants))) {
        if (!all(c("chr", "pos", "ref", "alt") %in% names(variants))) {
            stop("must supply either variant.id or chr,pos,ref,alt")
        }
    }
    x <- variantInfo(gds, alleles=TRUE, expanded=FALSE) %>%
        left_join(variants) %>%
        AnnotatedDataFrame()
    variantData(gds) <- x
    gds
}
