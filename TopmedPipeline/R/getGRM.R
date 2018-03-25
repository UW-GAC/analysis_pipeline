#' Load a GRM or kinship matrix
#'
#' @param config Config object (named vector) with params "pcrelate_file", "grm_file"
#' @param sample.id Vector of samples to include
#' @return GRM or kinship matrix. \code{NULL} if both pcrelate_file and grm_file are \code{NA}.
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
        pcr <- openfn.gds(config["pcrelate_file"])
        grm <- pcrelateMakeGRM(pcr, scan.include=sample.id, scaleKin=2)
        closefn.gds(pcr)
    } else if (!is.na(config["grm_file"])) {
        if (tools::file_ext(config["grm_file"]) == "gds") {
            x <- openfn.gds(config["grm_file"])
            samp <- read.gdsn(index.gdsn(x, "sample.id"))
            sel <- samp %in% sample.id
            grm <- readex.gdsn(index.gdsn(x, "grm"), sel=list(sel,sel))
            colnames(grm) <- rownames(grm) <- samp[sel]
            closefn.gds(x)
        } else {
            x <- getobj(config["grm_file"])
            if ("grm" %in% names(x)) {
                colnames(x$grm) <- rownames(x$grm) <- x$sample.id
                keep <- x$sample.id %in% sample.id
                grm <- x$grm[keep,keep]
            } else {
                keep <- colnames(x) %in% sample.id
                grm <- x[keep, keep]
            }
        }
    } else {
        grm <- NULL
    }
    return(grm)
}
