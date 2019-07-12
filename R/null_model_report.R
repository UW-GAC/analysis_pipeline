library(argparser)
library(TopmedPipeline)
sessionInfo()

argp <- arg_parser("Analysis report")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c()
optional <- c("out_file"="null_model_report")
config <- setConfigDefaults(config, required, optional)
print(config)

custom_render_markdown("null_model_report", config["out_file"])
