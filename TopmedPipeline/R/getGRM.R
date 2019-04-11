#' Load a GRM or kinship matrix
#'
#' Load a GRM or kinship matrix
#'
#' getGRM returns a Genetic Relationship Matrix from a pcrelate or grm file.
#'
#' @param config Config object (named vector) with params "pcrelate_file", "grm_file"
#' @param sample.id Vector of samples to include
#' @return List of GRMs or kinship matrices. \code{NULL} if all file names in config are \code{NA}.
#'
#' @export
getGRM <- function(config, sample.id=NULL) {
    if (!is.na(config["pcrelate_file"]) & !is.na(config["grm_file"])) {
        stop("Only one of pcrelate_file and grm_file may be specified")
    }
    
    ## load GRM for selected samples only
    if (!is.na(config["pcrelate_file"])) {
        files <- .splitFiles(config["pcrelate_file"])
        grm <- lapply(files, .readGRM, sample.id, matrix.name="kinship")
    } else if (!is.na(config["grm_file"])) {
        files <- .splitFiles(config["grm_file"])
        grm <- lapply(files, .readGRM, sample.id, matrix.name="grm")
    } else {
        grm <- NULL
    }

    return(grm)
}


.splitFiles <- function(f) {
    strsplit(f, " ", fixed=TRUE)[[1]]
}


#' @importFrom gdsfmt openfn.gds closefn.gds index.gdsn read.gdsn readex.gdsn
#' @import Matrix
#' @noRd
.readGRM <- function(f, sample.id, matrix.name="grm") {
    if (tools::file_ext(f) == "gds") {
        x <- openfn.gds(f)
        samp <- read.gdsn(index.gdsn(x, "sample.id"))
        if (is.null(sample.id)) sample.id <- samp
        sel <- samp %in% sample.id
        grm <- readex.gdsn(index.gdsn(x, matrix.name), sel=list(sel,sel))
        colnames(grm) <- rownames(grm) <- samp[sel]
        closefn.gds(x)
    } else {
        x <- getobj(f)
        if (matrix.name %in% names(x)) {
            colnames(x[[matrix.name]]) <- rownames(x[[matrix.name]]) <- x$sample.id
            if (!is.null(sample.id)) {
                keep <- x$sample.id %in% sample.id
                grm <- x[[matrix.name]][keep,keep]
            } else {
                grm <- x
            }
        } else {
            if (!is.null(sample.id)) {
                keep <- colnames(x) %in% sample.id
                grm <- x[keep, keep]
            } else {
                grm <- x
            }
        }
    }
    grm
}
