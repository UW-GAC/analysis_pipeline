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
    requireNamespace("SeqArray")
    stopifnot(length(gds.list) == length(id.list))
    n <- 0
    new.id <- list()
    for (i in seq_along(gds.list)) { 
        new.id[[i]] <- id.list[[i]] + n
        ni <- SeqArray::seqSummary(gds.list[[i]], "variant.id")
        n <- n + ni
    }
    stopifnot(all(id.list[[1]] == new.id[[1]]))
    seq.ids <- unlist(new.id)
    stopifnot(sum(duplicated(seq.ids)) == 0)
    return(seq.ids)
}

