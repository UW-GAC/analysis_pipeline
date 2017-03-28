#' Count number of threads available to the R process
#'
#' @return Number of threads
#'
#' @export
countThreads <- function() {
    nSlots <- Sys.getenv("NSLOTS")
    nThreads <- ifelse(is.na(strtoi(nSlots) >= 1), 1, strtoi(nSlots))
    if (nThreads == 0) nThreads <- 1
    message(paste("Running with", nThreads,"thread(s)."))
    nThreads
}
