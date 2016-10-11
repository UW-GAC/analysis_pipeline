library(argparser)
library(TopmedPipeline)
library(Biobase)
library(GENESIS)
library(gdsfmt)
sessionInfo()

argp <- arg_parser("Null model for association tests")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("outcome",
              "pca_file",
              "phenotype_file")
optional <- c("pcrelate_file"=NA,
              "binary"=FALSE,
              "covars"=NA,
              "inverse_normal"=FALSE,
              "n_pcs"=3,
              "out_file"="null_model.RData",
              "sample_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

# get phenotypes
annot <- getobj(config["phenotype_file"])

# select samples
if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
    annot <- annot[annot$sample.id %in% sample.id,]
} else {
    sample.id <- annot$sample.id
}

# get PCs
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

# outcome and covariates
outcome <- config["outcome"]
if (!is.na(config["covars"])) {
    covars <- strsplit(config["covars"], " ", fixed=TRUE)[[1]]
} else {
    covars <- NULL
}
covars <- c(covars, pccols)

if (as.logical(config["binary"])) {
    stopifnot(all(annot[[outcome]] %in% c(0,1)))
    family <- binomial
} else {
    family <- gaussian
}

if (!is.na(config["pcrelate_file"])) {
    ## load GRM for selected samples only
    pcr <- openfn.gds(config["pcrelate_file"])
    grm <- pcrelateMakeGRM(pcr, scan.include=sample.id, scaleKin=2)
    closefn.gds(pcr)

    message("Model: ", outcome, " ~ ", paste(c(covars, "(1|kinship)"), collapse=" + "))
    nullmod <- fitNullMM(annot, outcome=outcome, covars=covars,
                         covMatList=grm, scan.include=sample.id,
                         family=family)

    ## if we need an inverse normal transform, take residuals and refit null model
    ## for second model fit, use kinship but not other covariates
    if (as.logical(config["inverse_normal"])) {
        resid.norm <- rankNorm(nullmod$resid.marginal)
        annot$resid.norm <- resid.norm[match(annot$sample.id, nullmod$scanID)]
        message(paste0("resid.norm = rankNorm(resid.marginal(", outcome, " ~ ", paste(c(covars, "(1|kinship)"), collapse=" + "), "))"))
        ##message("Model: resid.norm ~ ", paste(c(pccols, "(1|kinship)"), collapse=" + "))
        message("Model: resid.norm ~ (1|kinship)")
        nullmod <- fitNullMM(annot, outcome="resid.norm", covars=NULL,
                             covMatList=grm, scan.include=sample.id,
                             family=family)
    }
} else {
    message("No kinship file specified, assuming samples are unrelated.")

    message("Model: ", outcome, " ~ ", paste(covars, collapse=" + "))
    nullmod <- fitNullReg(annot, outcome=outcome, covars=covars,
                         scan.include=sample.id, family=family)

    ## if we need an inverse normal transform, take residuals and refit null model
    ## for second model fit, use kinship but not other covariates
    if (as.logical(config["inverse_normal"])) {
        resid.norm <- rankNorm(nullmod$resid.response)
        annot$resid.norm <- resid.norm[match(annot$sample.id, nullmod$scanID)]
        message(paste0("resid.norm = rankNorm(resid.marginal(", outcome, " ~ ", paste(covars, collapse=" + "), "))"))
        message("Model: resid.norm")
        nullmod <- fitNullReg(annot, outcome="resid.norm", covars=NULL,
                              scan.include=sample.id, family=family)
    }
    
}

save(nullmod, file=config["out_file"])
