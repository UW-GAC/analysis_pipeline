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
optional <- c("out_prefix"="null_model",
			  "chromosomes"="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22")
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(basename(argv$config), ".null_model_fast_scoreSE.params"))

## need to be able to do this split by chr or not
## gds file can have two parts split by chromosome identifier
scorefile <- config["variant_score_file"]
bychrfile <- grepl(" ", scorefile) # do we have one file per chromosome?

if (bychrfile){
	chr <- strsplit(config["chromosomes"], " ", fixed=TRUE)[[1]]
	files <- sprintf(gsub(" ", "%s", config["variant_score_file"]), chr)
} else{
	files <- scorefile
}

# combine chr files
variantTable <- lapply(unname(files), getobj)

# get null model
nullModel <- getobj(config["null_model_file"])

# update the null model
nullModel <- nullModelFastScore(nullModel, variantTable)
outfile <- sprintf("%s.RData", config["out_prefix"])
save(nullModel, file = outfile)

# delete chr files
unlink(files)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
