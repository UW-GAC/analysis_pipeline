library(argparser)
library(TopmedPipeline)
sessionInfo()

argp <- arg_parser("Analysis report")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c()
optional <- c("out_file"="analysis_report")
config <- setConfigDefaults(config, required, optional)
print(config)

custom_render_markdown("analysis_report", config["out_file"])
