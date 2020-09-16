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
optional <- c("family"="gaussian",
              "inverse_normal"=TRUE,
              "out_prefix"="null_model",
              "n_categories_boxplot" = 10)
config <- setConfigDefaults(config, required, optional)
print(config)

print('rendering report for original null model')
p <- list(
  pipeline_version = argv$version,
  n_categories_boxplot = config["n_categories_boxplot"]
)
outfile <- sprintf("%s_report", config["out_prefix"])
custom_render_markdown("null_model_report", outfile, parameters = p)

# Generate the report for the inverse normal, if necessary.
if (as.logical(config["inverse_normal"]) & config["family"] == "gaussian") {
  print('rendering report for inverse normal transform')
  p <- c(p, list(invnorm = TRUE))
  outfile <- sprintf("%s_invnorm_report", config["out_prefix"])
  custom_render_markdown("null_model_report", outfile, parameters = p)
}
