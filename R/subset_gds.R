library(argparser)
library(TopmedPipeline)
library(SeqArray)
sessionInfo()

argp <- arg_parser("Subset GDS file")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("gds_file",
              "subset_gds_file")
optional <- c("sample_include_file"=NA,
              "variant_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

gds <- seqOpen(config["gds_file"])

if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
} else {
    sample.id <- NULL
}

if (!is.na(config["variant_include_file"])) {
    variant.id <- getobj(config["variant_include_file"])
} else {
    variant.id <- NULL
}

seqSetFilter(gds, sample.id=sample.id, variant.id=variant.id)
seqExport(gds, config["subset_gds_file"])

seqClose(gds)
