#' Get phenotypes
#'
#' @param config Config object (named vector) with params "phenotype_file", "n_pcs", "pca_file",
#'   "outcome", "covars", "sample_include_file", "gds_file", "conditional_variant_file"
#' @return list with annot, outcome, covars, sample.id
#'
#' @import Biobase
#' @importFrom stats complete.cases
#' @importFrom dplyr left_join
#' @export
getPhenotypes <- function(config) {

    ## get phenotypes
    annot <- getobj(config["phenotype_file"])

    ## get PCs
    n_pcs <- as.integer(config["n_pcs"])
    if (n_pcs > 0) {
        if (is.na(config["pca_file"])) stop("must specify pca_file if n_pcs > 0")
        pca <- getobj(config["pca_file"])
        if ("vectors" %in% names(pca)) {
            pca <- pca$vectors
        }
        pcs <- pca[,1:n_pcs,drop=FALSE]
        pccols <- paste0("PC", 1:n_pcs)
        colnames(pcs) <- pccols
        pcs <- pcs[match(annot$sample.id, rownames(pcs)),,drop=FALSE]
        rownames(pcs) <- NULL
        pData(annot) <- cbind(pData(annot), pcs)
    } else {
        pccols <- NULL
    }

    ## outcome and covariates
    outcome <- unname(config["outcome"])
    covars <- .parseParam(config["covars"])
    covars <- c(covars, pccols)
    group.var <- unname(config["group_var"])
    if (is.na(group.var)) group.var <- NULL

    annot <- annot[,unique(c("sample.id", outcome, covars, group.var))]

    ## select samples
    if (!is.na(config["sample_include_file"])) {
        sample.id <- getobj(config["sample_include_file"])
    } else {
        sample.id <- annot$sample.id
    }

    ## conditional variants
    if (!is.na(config["conditional_variant_file"])) {
        vars <- .conditionalVariants(config, sample.id=sample.id)
        pData(annot) <- left_join(pData(annot), vars, by="sample.id")
        covars <- c(covars, names(vars)[-1])
    }

    cc <- annot$sample.id[complete.cases(pData(annot))]
    sample.id <- intersect(sample.id, cc)

    ## match annot to gds file
    if (!is.na(config["gds_file"])) {
        annot <- .matchAnnotGdsConfig(config, annot)
    }

    list(annot=annot, outcome=outcome, covars=covars, group.var=group.var, sample.id=sample.id)
    
}


#' Parse space-separated parameter list
#'
#' @param param Parameter to parse
#' @return \code{NULL} if \code{param} is \code{NA}, vector otherwise
#'
#' @noRd
.parseParam <- function(param) {
    if (!is.na(param)) {
        strsplit(param, " ", fixed=TRUE)[[1]]
    } else {
        NULL
    }
}


.conditionalVariants <- function(config, sample.id) {
    dat <- getobj(config["conditional_variant_file"])
    if ("chromosome" %in% names(dat)) names(dat)[names(dat) == "chromosome"] <- "chr"
    stopifnot(all(c("chr", "variant.id") %in% names(dat)))
    geno <- do.call(cbind, lapply(unique(dat$chr), function(c) {
        vars <- dat$variant.id[dat$chr == c]
        gdsfile <- insertChromString(config["gds_file"], c)
        .genotypes(gdsfile, variant.id=vars, sample.id=sample.id)
    }))
    colnames(geno) <- paste0("var_", colnames(geno))
    data.frame(sample.id=rownames(geno), geno, row.names=1:nrow(geno), stringsAsFactors=FALSE)    
}


#' Return genotypes (alt dosage) for variants
#'
#' @param gdsfile Filename for GDS file
#' @param variant.id Vector of variant IDs
#' @param sample.id Vector of sample IDs
#' @return data.frame with genotypes for selected variants
#'
#' @import SeqArray
#' @importFrom SeqVarTools altDosage
#'
#' @noRd
.genotypes <- function(gdsfile, variant.id, sample.id=NULL) {
    gds <- seqOpen(gdsfile)
    seqSetFilter(gds, variant.id=variant.id, sample.id=sample.id, verbose=FALSE)
    geno <- altDosage(gds)
    seqClose(gds)
    geno
}


#' Return an expanded AnnotatedDataFrame with the same sample.id as the gds_file
#'
#' @param config  Config object (named vector) with params "gds_file"
#' @param annot AnnotatedDataFrame with sample.id matching gds file
#' @return AnnotatedDataFrame with sample.id matching gds file, and NA values for all
#'   samples not in original annot
#'
#' @import Biobase
#' @importFrom dplyr left_join
#'
#' @noRd
.matchAnnotGdsConfig <- function(config, annot) {
    tmp <- sub(" ", "[[:alnum:]]+", config["gds_file"])
    gdsfile <-  list.files(path=dirname(tmp), pattern=basename(tmp), full.names=TRUE)[1]
    gds <- seqOpen(gdsfile)
    annot <- matchAnnotGds(gds, annot)
    seqClose(gds)
    annot
}

#' Return an expanded AnnotatedDataFrame with the same sample.id as the gds_file
#'
#' @param gds SeqVarGDSClass object
#' @param annot AnnotatedDataFrame with sample.id matching gds file
#' @return AnnotatedDataFrame with sample.id matching gds file, and NA values for all
#'   samples not in original annot
#'
#' @import Biobase
#' @importFrom dplyr left_join
#'
#' @export
matchAnnotGds <- function(gds, annot) {
    sample.id <- seqGetData(gds, "sample.id")
    if (!(isTRUE(all.equal(sample.id, annot$sample.id)))) {
        dat <- data.frame(sample.id, stringsAsFactors=FALSE)
        pData(annot) <- left_join(dat, pData(annot), by="sample.id")
    }
    annot
}
