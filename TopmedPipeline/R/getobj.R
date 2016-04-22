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
