#' Return a small version of the null model (no NxN matrices)
#'
#' @param nullmod null model object
#' @return null model object without large matrices
#'
#' @export
smallNullModel <- function(nullmod) {
    nullmod$cholSigmaInv <- NULL
    nullmod$CX <- NULL
    nullmod$CXCXI <- NULL
    nullmod
}
