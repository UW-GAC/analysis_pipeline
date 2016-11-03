library(argparser)
library(TopmedPipeline)
sessionInfo()

argp <- arg_parser("LocusZoom plots")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("assoc_file",
              "assoc_type",
              "chromosome",
              "variant")
config <- setConfigDefaults(config, required, optional)
print(config)
