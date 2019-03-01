library(argparser)
library(TopmedPipeline)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Format KING results as Matrix")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("king_file")
optional <- c("kinship_threshold"=0.01104854, # 2^(-13/2), 5th degree
              "out_file"="king_Matrix.RData",
              "sample_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
} else {
    sample.id <- NULL
}

if (!is.na(config["kinship_threshold"])) {
    kin.thresh <- as.numeric(config["kinship_threshold"])
} else {
    kin.thresh <- NULL
}

mat <- kingToMatrix(king=config["king_file"],
                    sample.include=sample.id,
                    thresh=kin.thresh)

save(mat, file=config["out_file"])
