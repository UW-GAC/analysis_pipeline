#' Internal function to render markdown documents
#'
#' @param markdown_template_name The name of the markdown template, excluding ".Rmd"
#' @param output_file The name of the output file to create
#' @param parameters A named list of parameters to pass to the markdown document
#'
#' @seealso \code{\link[rmarkdown]{render}}
#'
#' @details
#' The markdown template is expected to live in the \code{inst/rmd} directory and an error
#' is raised if it is not found.
#'
#' @importFrom rmarkdown render
#' @export
custom_render_markdown <- function(markdown_template_name, output_file, parameters=NULL) {

    markdown_file <- file.path(system.file(package="TopmedPipeline"), "rmd",
                               sprintf("%s.Rmd", markdown_template_name))
    if (!file.exists(markdown_file)) {
        errmsg <- sprintf("markdown template not found: %s", markdown_file)
        stop(errmsg)
    }

    output_file <- sprintf("%s.Rmd", output_file)
    success <- file.copy(markdown_file, output_file, overwrite=TRUE)
    if (!success) {
        errmsg <- sprintf("could not copy markdown template to: %s", output_file)
        stop(errmsg)
    }

    rmarkdown::render(output_file, params=parameters, quiet=TRUE)
}

#' @name modelString
#' @rdname modelString
#' @title Create model strings for an analysis
#' @param outcome outcome string
#' @param covars vector of covariate strings
#' @param random vector of random effect strings
#' @param group_var grouping variable string
#' @param inverse_normal indicator of whether the model was refit with inverse-normal transformed residuals
#'
#' @return String representing analysis model
NULL

#' @describeIn modelString  Create a string for the entire model
#'
#' @export
modelString <- function(outcome, covars, random, group_var, inverse_normal = FALSE) {
    model.outcome <- modelOutcomeString(outcome, inverse_normal = inverse_normal)
    model.covars <- if (is.null(covars)) NULL else paste(covars, collapse=" + ")
    model.random <- if (is.null(random)) NULL else paste(paste0("(1|", random, ")"), collapse=" + ")
    model.var <- if (is.null(group_var)) NULL else paste0("var(", group_var, ")")
    model.string <- paste(c(model.covars, model.random, model.var), collapse=" + ")
    paste(model.outcome, model.string, sep=" ~ ")
}


#' @describeIn modelString Create model outcome strings
#'
#' @export
modelOutcomeString <- function(outcome, inverse_normal = FALSE) {
  if (inverse_normal) {
    return(sprintf("rankInvNorm(resid(%s))", outcome))
  } else {
    return(outcome)
  }
}
