#' Read a configuration file
#'
#' Functions for manipulating pipeline configuration files
#'
#' \code{readConfig} returns a named character vector of parameter values.
#' 
#' \code{writeConfig} writes a named character vector to a file.
#' 
#' \code{setConfigDefaults} takes a named character vector returned by
#' \code{readConfig} and adds additional parameters in \code{optional} with default values.
#' An error will result if a parameter in \code{required} is missing from \code{config}.
#'
#' @param file file where column 1 is parameter name and column 2 is value.
#' @param ... additional arguments to \code{read.table}
#' @examples
#' file <- tempfile()
#' write.table(cbind(letters[1:10], 1:10), file=file, quote=FALSE,
#'             row.names=FALSE, col.names=FALSE)
#' config <- readConfig(file)
#' 
#' required <- letters[1:5]
#' optional <- setNames(11:15, letters[11:15])
#' config <- setConfigDefaults(config, required, optional)
#' 
#' writeConfig(config, file)
#' 
#' unlink(file)
#'
#' @importFrom utils read.table
#' @export
readConfig <- function(file, ...) {
  config.table <- read.table(file, as.is=TRUE, ...)
  if (any(duplicated(config.table[, 1]))) stop("duplicated parameters in config file are not allowed!")
  config <- config.table[,2]
  names(config) <- config.table[,1]
  # recode tabs
  config[config %in% "\\t"] <- "\t"
 
  return(config)
}

#' @param config named character vector
#' @rdname readConfig
#'
#' @importFrom utils write.table
#' @export
writeConfig <- function(config, file, ...) {
  write.table(config, file=file, col.names=FALSE, ...)
}

#' @param required character vector of required parameter names
#' @param optional named vector of optional parameter values
#' @rdname readConfig
#'
#' @export
setConfigDefaults <- function(config, required, optional) {
  # optional is a named list of default values
  default <- unname(optional)
  optional <- names(optional)
  
  config.params <- names(config)
  found.params <- intersect(config.params, c(required, optional))
  if (length(found.params) > 0) {
    message("found parameters: ", paste(found.params, collapse=", "))
  }
  
  # if required params not in config, stop
  missing.params <- setdiff(required, config.params)
  if (length(missing.params) > 0) {
    stop("missing required parameters: ", paste(missing.params, collapse=", "))
  }

  # if not in config, set default value
  set.params <- setdiff(optional, config.params)
  if (length(set.params) > 0) {
    config[set.params] <- default[match(set.params, optional)]
    message("using default values: ", paste(set.params, collapse=", "))
  }
  
  # note unsed params in config
  extra.params <- setdiff(config.params, c(required, optional))
  if (length(extra.params) > 0) {
    message("unused parameters: ", paste(extra.params, collapse=", "))
  }
  
  # return config with default values set
  config <- config[c(required, optional)]
  return(config)
}
