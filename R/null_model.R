library(argparser)
library(TopmedPipeline)
library(Biobase)
library(genesis2)
library(genesis2tests)
library(gdsfmt)
sessionInfo()

argp <- arg_parser("Null model for association tests")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("outcome",
              "phenotype_file")
optional <- c("gds_file"=NA, # required for conditional variants
              "pca_file"=NA,
              "pcrelate_file"=NA,
              "grm_file"=NA,
              "binary"=FALSE,
              "conditional_variant_file"=NA,
              "covars"=NA,
              "group_var"=NA,
              "inverse_normal"=TRUE,
              "n_pcs"=3,
              "out_file"="null_model.RData",
              "out_phenotype_file"="phenotypes.RData",
              "rescale_variance"="marginal",
              #"resid_covars"=TRUE,
              "sample_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(basename(argv$config), ".null_model.params"))

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
nullmod <- fitNullModel2(annot, outcome=outcome, covars=covars,
                         cov.mat=grm, sample.id=sample.id,
                         family=family, group.var=group.var)

## if we need an inverse normal transform, take residuals and refit null model
if (as.logical(config["inverse_normal"])) {
    if (is.null(group.var)) {
        rankNorm.option <- "all"
        rescale <- "None"
    } else {
        rankNorm.option <- "by.group"        
        if (config["rescale_variance"] == "varcomp") {
            rescale <- "model"
        } else if (config["rescale_variance"] == "marginal") {
            rescale <- "residSD"
        }
    }
    
    nullmod <- updateNullModOutcome(nullmod,
                                    covMatList=grm,
                                    rankNorm.option=rankNorm.option,
                                    rescale=rescale)
}

save(nullmod, file=config["out_file"])

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
