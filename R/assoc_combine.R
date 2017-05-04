library(argparser)
library(TopmedPipeline)
sessionInfo()

argp <- arg_parser("Association test - combine files")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argv <- parse_args(argp)
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)

required <- c("assoc_type",
              "out_prefix")
optional <- c()
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(argv$config, ".assoc_combine.params"))

file.pattern <- paste0(basename(config["out_prefix"]), "_chr", chr, "_seg[[:digit:]]+.RData")
files <- list.files(path=dirname(config["out_prefix"]), pattern=file.pattern, full.names=TRUE)

assoc <- combineAssoc(files, config["assoc_type"])

save(assoc, file=constructFilename(config["out_prefix"], chr))

# delete segment files
unlink(files)
