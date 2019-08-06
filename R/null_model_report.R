library(argparser)
library(TopmedPipeline)
sessionInfo()

argp <- arg_parser("Analysis report")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("out_prefix")
optional <- c("binary"=FALSE,
              "inverse_normal"=TRUE)
config <- setConfigDefaults(config, required, optional)
print(config)

print('rendering report for original null model')
p <- list(pipeline_version = argv$version)
outfile <- sprintf("%s_report", config["out_prefix"])
custom_render_markdown("null_model_report", outfile, parameters = p)

# Generate the report for the inverse normal, if necessary.
if (as.logical(config["inverse_normal"]) & !as.logical(config["binary"])) {
  print('rendering report for inverse normal transform')
  p <- c(p, list(invnorm = TRUE))
  outfile <- sprintf("%s_invnorm_report", config["out_prefix"])
  custom_render_markdown("null_model_report", outfile, parameters = p)
}
