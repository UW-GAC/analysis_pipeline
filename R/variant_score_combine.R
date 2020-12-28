library(argparser)
library(TopmedPipeline)
sessionInfo()

argp <- arg_parser("Variant score and SE - combine files")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("out_prefix")
optional <- c()
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(basename(argv$config), ".variant_score_combine.params"))

file.pattern <- paste0(basename(config["out_prefix"]), "_chr[[:digit:]]+.RData")
files <- list.files(path=dirname(config["out_prefix"]), pattern=file.pattern, full.names=TRUE)
chrs <- sub("_chr", "", regmatches(basename(files), regexpr("_chr[[:digit:]]+", basename(files))))
files <- files[order(as.integer(chrs))]

# combine chr files
tab <- data.table::rbindlist(lapply(unname(files), getobj))

outfile <- sprintf("%s.RData", config["out_prefix"])
save(tab, file = outfile)

# delete chr files
unlink(files)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
