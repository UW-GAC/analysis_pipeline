#' Load a GRM or kinship matrix
#'
#' Load a GRM or kinship matrix
#'
#' getGRM returns a Genetic Relationship Matrix from a pcrelate or grm file.
#'
#' getKinship returns a kinship matrix from a pcrelate or king file.
#'
#' @param config Config object (named vector) with params "pcrelate_file", "grm_file" or "king_file"
#' @param sample.id Vector of samples to include
#' @return List of GRMs or kinship matrices. \code{NULL} if all file names in config are \code{NA}.
#'
#' @export
getGRM <- function(config, sample.id) {
    if (!is.na(config["pcrelate_file"]) & !is.na(config["grm_file"])) {
        stop("Only one of pcrelate_file and grm_file may be specified")
    }
    
    if (!is.na(config["pcrelate_file"])) {
        ## load GRM for selected samples only
        files <- .splitFiles(config["pcrelate_file"])
        grm <- lapply(files, .readPCR, sample.id, scaleKin=2)
    } else if (!is.na(config["grm_file"])) {
        files <- .splitFiles(config["grm_file"])
        grm <- lapply(files, .readGRM, sample.id, matrix.name="grm")
    } else {
        grm <- NULL
    }

    return(grm)
}

#' @rdname getGRM
#'
#' @export
getKinship <- function(config, sample.id) {
    if (!is.na(config["pcrelate_file"]) & !is.na(config["king_file"])) {
        stop("Only one of pcrelate_file and king_file may be specified")
    }
    
    if (!is.na(config["pcrelate_file"])) {
        ## load GRM for selected samples only
        files <- .splitFiles(config["pcrelate_file"])
        grm <- lapply(files, .readPCR, sample.id, scaleKin=1)
    } else if (!is.na(config["king_file"])) {
        files <- .splitFiles(config["king_file"])
        grm <- lapply(files, .readGRM, sample.id, matrix.name="kinship")
    } else {
        grm <- NULL
    }

    return(grm)
}

.splitFiles <- function(f) {
    strsplit(f, " ", fixed=TRUE)[[1]]
}

#' @importFrom GENESIS pcrelateMakeGRM
#' @noRd
.readPCR <- function(f, sample.id, scaleKin=2) {
    if (tools::file_ext(f) == "gds") {
        pcr <- openfn.gds(f)
        grm <- pcrelateMakeGRM(pcr, scan.include=sample.id, scaleKin=scaleKin)
        closefn.gds(pcr)
    } else {
        pcr <- getobj(f)
        grm <- pcrelateMakeGRM(pcr, scan.include=sample.id, scaleKin=scaleKin)
    }
    grm
}

#' @importFrom gdsfmt openfn.gds closefn.gds index.gdsn read.gdsn readex.gdsn
#' @import Matrix
#' @noRd
.readGRM <- function(f, sample.id, matrix.name="grm") {
    if (tools::file_ext(f) == "gds") {
        x <- openfn.gds(f)
        samp <- read.gdsn(index.gdsn(x, "sample.id"))
        sel <- samp %in% sample.id
        grm <- readex.gdsn(index.gdsn(x, matrix.name), sel=list(sel,sel))
        colnames(grm) <- rownames(grm) <- samp[sel]
        closefn.gds(x)
    } else {
        x <- getobj(f)
        if (matrix.name %in% names(x)) {
            colnames(x[[matrix.name]]) <- rownames(x[[matrix.name]]) <- x$sample.id
            keep <- x$sample.id %in% sample.id
            grm <- x[[matrix.name]][keep,keep]
        } else {
            keep <- colnames(x) %in% sample.id
            grm <- x[keep, keep]
        }
    }
    grm
}
