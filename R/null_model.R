library(argparser)
library(TopmedPipeline)
library(Biobase)
library(GENESIS)
library(gdsfmt)
sessionInfo()

argp <- arg_parser("Null model for association tests")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("outcome",
              "phenotype_file",
              "out_prefix")
optional <- c("gds_file"=NA, # required for conditional variants
              "pca_file"=NA,
              "pcrelate_file"=NA,
              "grm_file"=NA,
              "binary"=FALSE,
              "conditional_variant_file"=NA,
              "covars"=NA,
              "group_var"=NA,
              "inverse_normal"=TRUE,
              "n_pcs"=0,
              "norm_bygroup"=FALSE,
              "out_phenotype_file"="phenotypes.RData",
              "rescale_variance"="marginal",
              "sample_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(basename(argv$config), ".null_model.params"))

# get the number of threads available
# this should speed up matrix calculations if we are running parallel MKL
countThreads()

# get phenotypes
phen <- getPhenotypes(config)
annot <- phen[["annot"]]
outcome <- phen[["outcome"]]
covars <- phen[["covars"]]
group.var <- phen[["group.var"]]
sample.id <- phen[["sample.id"]]

save(annot, file=config["out_phenotype_file"])

if (as.logical(config["binary"])) {
    stopifnot(all(annot[[outcome]] %in% c(0,1,NA)))
    family <- binomial
} else {
    family <- gaussian
}

# kinship matrix or GRM
grm <- getGRM(config, sample.id)

# print model
random <- if (!is.na(config["pcrelate_file"])) "kinship" else if (!is.na(config["grm_file"])) "GRM" else NULL
model.string <- modelString(outcome, covars, random, group.var)
message("Model: ", model.string)
message(length(sample.id), " samples")

## fit null model allowing heterogeneous variances among studies
nullmod <- fitNullModel(annot, outcome=outcome, covars=covars,
                        cov.mat=grm, sample.id=sample.id,
                        family=family, group.var=group.var)

# Save a smaller version of the original null model.
nullmod_small <- smallNullModel(nullmod)
outfile <- sprintf("%s_small.RData", config["out_prefix"])
save(nullmod_small, file = outfile)


## if we need an inverse normal transform, take residuals and refit null model
if (as.logical(config["inverse_normal"]) & !as.logical(config["binary"])) {

    if (as.logical(config["norm_bygroup"]) & !is.null(group.var)) {
        norm.option <- "by.group"
    } else {
        norm.option <- "all"
    }
    if (config["rescale_variance"] == "varcomp") {
        rescale <- "model"
    } else if (config["rescale_variance"] == "marginal") {
        rescale <- "residSD"
    } else {
        rescale <- "none"
    }

    nullmod <- nullModelInvNorm(nullmod, cov.mat=grm,
                                norm.option=norm.option,
                                rescale=rescale)

    # Save a smaller version of the null model.
    nullmod_small <- smallNullModel(nullmod)
    outfile <- sprintf("%s_invnorm_small.RData", config["out_prefix"])
    save(nullmod_small, file = outfile)

    # change filename to indicate invnorm
    outfile <- sprintf("%s_invnorm.RData", config["out_prefix"])
}

# save full version of final model
save(nullmod, file = outfile)



# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
