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
              "phenotype_file")
optional <- c("gds_file"=NA, # required for conditional variants
              "pca_file"=NA,
              "relatedness_matrix_file"=NA,
              "family"="gaussian",
              "conditional_variant_file"=NA,
              "covars"=NA,
              "group_var"=NA,
              "inverse_normal"=TRUE,
              "n_pcs"=0,
              "norm_bygroup"=FALSE,
              "out_phenotype_file"="phenotypes.RData",
              "out_prefix"="null_model",
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

family <- config["family"]
stopifnot(family %in% c("gaussian", "binomial", "poisson"))
if (family == "binomial") {
    stopifnot(all(annot[[outcome]] %in% c(0,1,NA)))
} else if (family == "poisson") {
    stopifnot(any(na.omit(annot[[outcome]]) < 0))
}

# kinship matrix or GRM
grm <- getGRM(config, sample.id)

# print model
random <- if (!is.na(config["relatedness_matrix_file"])) "relatedness" else NULL
model.string <- modelString(outcome, covars, random, group.var)
message("Model: ", model.string)
message(length(sample.id), " samples")

## fit null model allowing heterogeneous variances among studies
nullmod <- fitNullModel(annot, outcome=outcome, covars=covars,
                        cov.mat=grm, sample.id=sample.id,
                        family=family, group.var=group.var)

# Add the model string as a temporary fix until it can be added to GENESIS null models.
nullmod$model.string <- model.string

# Save a smaller version of the original null model.
nullmod_small <- smallNullModel(nullmod)
outfile <- sprintf("%s_reportonly.RData", config["out_prefix"])
save(nullmod_small, file = outfile)

# outfile if there is no invnorm
outfile <- sprintf("%s.RData", config["out_prefix"])


## if we need an inverse normal transform, take residuals and refit null model
if (as.logical(config["inverse_normal"]) & family == "gaussian") {

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

    model.string <- modelString(outcome, covars, random, group.var,
                                inverse_normal = TRUE)
    message("Refitting model: ", model.string)
    message(length(sample.id), " samples")

    nullmod <- nullModelInvNorm(nullmod, cov.mat=grm,
                                norm.option=norm.option,
                                rescale=rescale)

    # Update the model string so it has the inverse normal outcome.
    # This is a temporary fix until the model.string can be added to GENESIS null models.
    nullmod$model.string <- model.string

    # Save a smaller version of the null model.
    nullmod_small <- smallNullModel(nullmod)
    outfile <- sprintf("%s_invnorm_reportonly.RData", config["out_prefix"])
    save(nullmod_small, file = outfile)

    # change filename to indicate invnorm
    outfile <- sprintf("%s_invnorm.RData", config["out_prefix"])
}

# save full version of final model
save(nullmod, file = outfile)



# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
