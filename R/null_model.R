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
              "group_var"=NA,
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

# heterogeneous residual variance
if (!is.na(config["group_var"])) {
    group.var <- config["group_var"]
} else {
    group.var <- NULL
}

if (!is.na(config["pcrelate_file"])) {
    ## load GRM for selected samples only
    pcr <- openfn.gds(config["pcrelate_file"])
    grm <- pcrelateMakeGRM(pcr, scan.include=sample.id, scaleKin=2)
    closefn.gds(pcr)

    
    ## fit null model allowing heterogeneous variances among studies
    message("Model: ", outcome, " ~ ", paste(c(covars, "(1|kinship)"), collapse=" + "))
    nullmod <- fitNullMM(annot, outcome=outcome, covars=covars,
                         covMatList=grm, scan.include=sample.id,
                         family=family, group.var=group.var)

    ## if we need an inverse normal transform, take residuals and refit null model
    ## for second model fit, use kinship but not other covariates
    if (as.logical(config["inverse_normal"])) {
        if (is.null(group.var)) {
            annot <- addInvNorm(annot, nullmod, outcome, covars)
        } else {
            groups <- unique(annot[[group.var]])
            ## inverse-normal transform residuals from each study separately (mean=0, var=1)
            ## rescale the inverse-normal residuals to have study-specific variances = kinship variance component + study-specific residual
            resid.group <- do.call(rbind, lapply(groups, function(g) {
                samp.g <- intersect(nullmod$scanID, annot$sample.id[annot[[group.var]] == g])
                resid.g <- nullmod$resid.marginal[nullmod$scanID %in% samp.g]
                resid.scale <- nullmod$varComp["V_A"] + nullmod$varComp[paste0("V_", g)]
                resid.norm <- rankNorm(resid.g) * sqrt(resid.scale)
                data.frame(sample.id=samp.g, resid.norm, stringsAsFactors=FALSE)
            }))
            annot$resid.norm <- resid.group$resid.norm[match(annot$sample.id, resid.group$sample.id)]
        }
        
        ## fit null model again with these residuals as outcome and allowing heterogeneous variances
        nullmod <- fitNullMM(annot, outcome="resid.norm", covars=NULL,
                             covMatList=grm, scan.include=sample.id,
                             family=family, group.var=group.var)
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
