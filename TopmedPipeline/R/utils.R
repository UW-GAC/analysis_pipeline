#load an Rdata file and return the object
#(using load directly returns only the object name)
getobj <- function(Rdata) {
  objname <- load(Rdata)
  if (length(objname) > 1) {
    warning(paste("Multiple objects stored in file", Rdata,
                  "\nReturning only the first object"))
  }
  return(get(objname))
}


insertChromString <- function(x, chr, err=NULL) {
    if (!is.null(err) & !(grepl(" ", x, fixed=TRUE))) {
        stop(paste(err, "must have a blank space to insert chromosome number"))
    }
    sub(" ", chr, x, fixed=TRUE)
}


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


calculateLambda <- function(stat, df) {
    if (any(sum(stat < 0, na.rm=TRUE)))
        stop("no negative values allowed in stat (does beta/se need to be squared?)")
    median(stat, na.rm=TRUE) / qchisq(0.5, df=df)
}


rankNorm <- function(x) qnorm((rank(x) - 0.5)/length(x))


intToChr <- function(chr) {
    if (is.na(chr)) return(NA)
    if (chr == 23) return("X")
    if (chr == 24) return("Y")
    as.character(chr)
}


constructFilename <- function(prefix, chromosome=NA, segment=NA) {
    chr <- if (is.na(chromosome)) "" else paste0("_chr", chromosome)
    seg <- if (is.na(segment)) "" else paste0("_seg", segment)
    paste0(prefix, chr, seg, ".RData")
}
