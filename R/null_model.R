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
phen <- getPhenotypes(config)
annot <- phen[["annot"]]
outcome <- phen[["outcome"]]
covars <- phen[["covars"]]
sample.id <- phen[["sample.id"]]

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
        annot <- addInvNorm(annot, nullmod, outcome, covars)
        nullmod <- fitNullMM(annot, outcome="resid.norm", covars=NULL,
                             covMatList=grm, scan.include=sample.id,
                             family=family)
    }
} else {
    message("No kinship file specified, assuming samples are unrelated.")

    message("Model: ", outcome, " ~ ", paste(covars, collapse=" + "))
    nullmod <- fitNullReg(annot, outcome=outcome, covars=covars,
                         scan.include=sample.id, family=family)

    if (as.logical(config["inverse_normal"])) {
        annot <- addInvNorm(annot, nullmod, outcome, covars)
        nullmod <- fitNullReg(annot, outcome="resid.norm", covars=NULL,
                              scan.include=sample.id, family=family)
    }
    
}

save(nullmod, file=config["out_file"])
