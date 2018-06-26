#' Get an R object stored in an Rdata file
#'
#' Returns an R object stored in an Rdata file
#'
#' Loads an R object and stores it under a new name without creating a duplicate copy.
#' If multiple objects are stored in the same file, only the first one will be returned.
#' 
#' @param Rdata path to an Rdata file containing a single R object to load
#' @return The R object stored in \code{Rdata}
#' @examples
#' x <- 1:10
#' file <- tempfile()
#' save(x, file=file)
#' y <- getobj(file)
#' unlink(file)
#'
#' @export
getobj <- function(Rdata) {
  objname <- load(Rdata)
  if (length(objname) > 1) {
    warning(paste("Multiple objects stored in file", Rdata,
                  "\nReturning only the first object"))
  }
  return(get(objname))
}


#' Format a string by inserting chromosome into a blank space
#'
#' @param x Character string
#' @param chr Chromosome number (or character) to instert
#' @param err If not \code{NULL}, print this string with an error message about requiring a blank space
#' @return String \code{x} with \code{chr} inserted into blank space
#'
#' @export
insertChromString <- function(x, chr, err=NULL) {
    if (!is.null(err) & !(grepl(" ", x, fixed=TRUE))) {
        stop(paste(err, "must have a blank space to insert chromosome number"))
    }
    sub(" ", chr, x, fixed=TRUE)
}


#' Generate sequential variant ids
#'
#' Generate sequential variant ids
#'
#' For a set of GDS files each with ids numbered 1:n, convert to sequential ids
#' with 1:N over the combined set of files.
#' 
#' @param gds.list List of \code{\link[SeqArray]{SeqVarGDSClass}} objects
#' @param id.list List of vectors of variant ids corresponding to \code{gds.list}
#' @return Vector of sequential variant ids
#'
#' @import SeqArray
#' @export
sequentialVariantIds <- function(gds.list, id.list) {
    stopifnot(length(gds.list) == length(id.list))
    n <- 0
    new.id <- list()
    for (i in seq_along(gds.list)) { 
        new.id[[i]] <- id.list[[i]] + n
        ni <- seqSummary(gds.list[[i]], "variant.id")
        n <- n + ni
    }
    stopifnot(all(id.list[[1]] == new.id[[1]]))
    seq.ids <- unlist(new.id)
    stopifnot(sum(duplicated(seq.ids)) == 0)
    return(seq.ids)
}


#' Calculate lambda, the genomic inflation factor
#'
#' @param stat Vector of test statistics
#' @param df Degrees of freedom
#' @return lambda
#'
#' @importFrom stats median qchisq
#' @export
calculateLambda <- function(stat, df) {
    if (any(sum(stat < 0, na.rm=TRUE)))
        stop("no negative values allowed in stat (does beta/se need to be squared?)")
    median(stat, na.rm=TRUE) / qchisq(0.5, df=df)
}


#' Integer chromosome code to character
#'
#' @param chr Integer chromosome code
#' @return Character chromosome (23 -> X, 24 -> Y)
#'
#' @export
intToChr <- function(chr) {
    if (is.na(chr)) return(NA)
    if (chr == 23) return("X")
    if (chr == 24) return("Y")
    as.character(chr)
}


#' Construct filename from chromosome and segment
#'
#' @param prefix Prefix of filename
#' @param chromosome chromosome
#' @param segment segment
#' @return Character string of form "<prefix>_chr<chromosome>_seg<segment>.RData"
#'
#' @export
constructFilename <- function(prefix, chromosome=NA, segment=NA) {
    chr <- if (is.na(chromosome)) "" else paste0("_chr", chromosome)
    seg <- if (is.na(segment)) "" else paste0("_seg", segment)
    paste0(prefix, chr, seg, ".RData")
}
