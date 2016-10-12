getPhenotypes <- function(config) {

    ## get phenotypes
    annot <- getobj(config["phenotype_file"])

    ## select samples
    if (!is.na(config["sample_include_file"])) {
        sample.id <- getobj(config["sample_include_file"])
        annot <- annot[annot$sample.id %in% sample.id,]
    } else {
        sample.id <- annot$sample.id
    }

    ## get PCs
    n_pcs <- as.integer(config["n_pcs"])
    if (n_pcs > 0) {
        pca <- getobj(config["pca_file"])
        pcs <- pca$vectors[,1:n_pcs,drop=FALSE]
        pccols <- paste0("PC", 1:n_pcs)
        colnames(pcs) <- pccols
        sample.id <- intersect(sample.id, rownames(pcs))
        annot <- annot[annot$sample.id %in% sample.id,]
        pData(annot) <- cbind(pData(annot), pcs[as.character(sample.id),,drop=FALSE])
    } else {
        pccols <- NULL
    }

    ## outcome and covariates
    outcome <- config["outcome"]
    if (!is.na(config["covars"])) {
        covars <- strsplit(config["covars"], " ", fixed=TRUE)[[1]]
    } else {
        covars <- NULL
    }
    covars <- c(covars, pccols)

    annot <- annot[,c("sample.id", outcome, covars)]

    list(annot=annot, outcome=outcome, covars=covars)
    
}
