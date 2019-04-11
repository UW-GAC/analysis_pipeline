library(argparser)
library(TopmedPipeline)
sessionInfo()

argp <- arg_parser("duplicate discordance")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("out_prefix")
optional <- c()
config <- setConfigDefaults(config, required, optional)
print(config)

file.pattern <- paste0(basename(config["out_prefix"]), "_seg[[:digit:]]+.RData")
files <- list.files(path=dirname(config["out_prefix"]), pattern=file.pattern, full.names=TRUE)

res <- do.call(rbind, lapply(files, getobj))

save(res, file=constructFilename(config["out_prefix"]))

# delete segment files
unlink(files)
