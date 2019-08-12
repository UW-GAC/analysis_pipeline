library(argparser)
library(TopmedPipeline)
sessionInfo()

argp <- arg_parser("Association test - combine files")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)

required <- c("assoc_type",
              "out_prefix")
optional <- c("conditional_variant_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(basename(argv$config), ".assoc_combine.params"))

file.pattern <- paste0(basename(config["out_prefix"]), "_chr", chr, "_seg[[:digit:]]+.RData")
files <- list.files(path=dirname(config["out_prefix"]), pattern=file.pattern, full.names=TRUE)
segments <- sub("_seg", "", regmatches(basename(files), regexpr("_seg[[:digit:]]+", basename(files))))
files <- files[order(as.integer(segments))]

assoc <- combineAssoc(files, config["assoc_type"])

if (!is.na(config["conditional_variant_file"])) {
    assoc <- removeConditional(assoc, config["conditional_variant_file"])
}

save(assoc, file=constructFilename(config["out_prefix"], chr))

# delete segment files
unlink(files)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
