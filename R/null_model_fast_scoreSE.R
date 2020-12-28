library(argparser)
library(TopmedPipeline)
library(Biobase)
library(SeqVarTools)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Update null model for fast score approximation association tests")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("null_model_file",
			  "variant_score_file")
optional <- c("out_prefix"="null_model")
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(basename(argv$config), ".null_model_fast_scoreSE.params"))

# get null model
nullModel <- getobj(config["null_model_file"])

# get variant score table
variantTable <- getobj(config["variant_score_file"])

# update the null model
nullModel <- nullModelFastScore(nullModel, variantTable)
outfile <- sprintf("%s_fast_scoreSE.RData", config["out_prefix"])
save(nullModel, file = outfile)


# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
