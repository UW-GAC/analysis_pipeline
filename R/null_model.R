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
optional <- c("gds_file"=NA, # required for conditional variants
              "pcrelate_file"=NA,
              "grm_file"=NA,
              "binary"=FALSE,
              "conditional_chrom"=NA,
              "conditional_variants"=NA,
              "covars"=NA,
              "group_var"=NA,
              "inverse_normal"=TRUE,
              "n_pcs"=3,
              "out_file"="null_model.RData",
              "rescale_variance"=TRUE,
              "resid_covars"=TRUE,
              "sample_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(argv$config, ".null_model.params"))

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

# kinship matrix or GRM
grm <- getGRM(config, sample.id)

if (!is.null(grm)) {
    
    ## fit null model allowing heterogeneous variances among studies
    message("Model: ", outcome, " ~ ", paste(c(covars, "(1|kinship)"), collapse=" + "))
    nullmod <- fitNullMM(annot, outcome=outcome, covars=covars,
                         covMatList=grm, scan.include=sample.id,
                         family=family, group.var=group.var)

    ## if we need an inverse normal transform, take residuals and refit null model
    if (as.logical(config["inverse_normal"])) {
        if (is.null(group.var)) {
            annot <- addInvNorm(annot, nullmod, outcome, covars)
        } else {
            groups <- unique(annot[[group.var]])
            ## inverse-normal transform residuals from each study separately (mean=0, var=1)
            resid.group <- do.call(rbind, lapply(groups, function(g) {
                samp.g <- intersect(nullmod$scanID, annot$sample.id[annot[[group.var]] == g])
                resid.g <- nullmod$resid.marginal[nullmod$scanID %in% samp.g]
                resid.norm <- rankNorm(resid.g)
                ## rescale the inverse-normal residuals to have study-specific variances =
                ## kinship variance component + study-specific residual
                if (as.logical(config["rescale_variance"])) {
                    resid.scale <- nullmod$varComp["V_A"] + nullmod$varComp[paste0("V_", g)]
                    resid.norm <- resid.norm * sqrt(resid.scale)
                }
                data.frame(sample.id=samp.g, resid.norm, stringsAsFactors=FALSE)
            }))
            annot$resid.norm <- resid.group$resid.norm[match(annot$sample.id, resid.group$sample.id)]
        }
        
        ## fit null model again with these residuals as outcome and allowing heterogeneous variances
        resid.covars <- if (config["resid_covars"]) covars else NULL
        nullmod <- fitNullMM(annot, outcome="resid.norm", covars=resid.covars,
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
        resid.covars <- if (config["resid_covars"]) covars else NULL
        nullmod <- fitNullReg(annot, outcome="resid.norm", covars=resid.covars,
                              scan.include=sample.id, family=family)
    }
    
}

save(nullmod, file=config["out_file"])
