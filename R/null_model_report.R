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
optional <- c("out_file"="null_model_report",
              "binary"=FALSE,
              "inverse_normal"=TRUE)
config <- setConfigDefaults(config, required, optional)
print(config)

print('rendering report for original null model')
custom_render_markdown("null_model_report", config["out_file"])

# Generate the report for the inverse normal, if necessary.
if (as.logical(config["inverse_normal"]) & !as.logical(config["binary"])) {
  print('rendering report for inverse normal transform')
  outfile <- gsub("report", "invnorm_report", config["out_file"])
  print(sprintf("output file: %s", outfile))
  p <- list(invnorm = TRUE)
  custom_render_markdown("null_model_report", outfile, parameters = p)
}
