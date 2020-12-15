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
