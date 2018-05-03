#' Load a GRM or kinship matrix
#'
#' @param config Config object (named vector) with params "pcrelate_file", "grm_file"
#' @param sample.id Vector of samples to include
#' @return List of GRMs or kinship matrices. \code{NULL} if both pcrelate_file and grm_file are \code{NA}.
#'
#' @importFrom GENESIS pcrelateMakeGRM
#' @importFrom gdsfmt openfn.gds closefn.gds index.gdsn read.gdsn readex.gdsn
#' @import Matrix
#' @export
getGRM <- function(config, sample.id) {
    if (!is.na(config["pcrelate_file"]) & !is.na(config["grm_file"])) {
        stop("Only one of pcrelate_file and grm_file may be specified")
    }
    
    if (!is.na(config["pcrelate_file"])) {
        ## load GRM for selected samples only
        files <- .splitFiles(config["pcrelate_file"])
        grm <- lapply(files, .readPCR, sample.id)
    } else if (!is.na(config["grm_file"])) {
        files <- .splitFiles(config["grm_file"])
        grm <- lapply(files, .readGRM, sample.id)
    } else {
        grm <- NULL
    }

    return(grm)
}

.splitFiles <- function(f) {
    strsplit(f, " ", fixed=TRUE)[[1]]
}

.readPCR <- function(f, sample.id) {
    pcr <- openfn.gds(f)
    grm <- pcrelateMakeGRM(pcr, scan.include=sample.id, scaleKin=2)
    closefn.gds(pcr)
    grm
}

.readGRM <- function(f, sample.id) {
    if (tools::file_ext(f) == "gds") {
        x <- openfn.gds(f)
        samp <- read.gdsn(index.gdsn(x, "sample.id"))
        sel <- samp %in% sample.id
        grm <- readex.gdsn(index.gdsn(x, "grm"), sel=list(sel,sel))
        colnames(grm) <- rownames(grm) <- samp[sel]
        closefn.gds(x)
    } else {
        x <- getobj(f)
        if ("grm" %in% names(x)) {
            colnames(x$grm) <- rownames(x$grm) <- x$sample.id
            keep <- x$sample.id %in% sample.id
            grm <- x$grm[keep,keep]
        } else {
            keep <- colnames(x) %in% sample.id
            grm <- x[keep, keep]
        }
    }
    grm
}
