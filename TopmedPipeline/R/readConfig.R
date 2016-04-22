readConfig <- function(file, ...) {
  config.table <- read.table(file, as.is=TRUE, ...)
  if (any(duplicated(config.table[, 1]))) stop("duplicated parameters in config file are not allowed!")
  config <- config.table[,2]
  names(config) <- config.table[,1]
  # recode tabs
  config[config %in% "\\t"] <- "\t"
 
  return(config)
}

writeConfig <- function(config, file, ...) {
  write.table(config, file=file, col.names=FALSE, ...)
}

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
  return(config)
}
