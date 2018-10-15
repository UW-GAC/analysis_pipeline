library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Partition samples into related and unrelated sets")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("king_file")
optional <- c("kinship_file"=NA,
              "kinship_method"="king",
              "kinship_threshold"=0.04419417, # 2^(-9/2), 3rd degree
              "out_related_file"="related.RData",
              "out_unrelated_file"="unrelated.RData",
              #"pcrelate_file"=NA,
              "sample_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
    message("Using ", length(sample.id), " samples")
} else {
    sample.id <- NULL
    message("Using all samples")
}

# always use king estimates for ancestry divergence
# getKinship returns a list
divMat <- getKinship(config["king_file"], sample.id)[[1]]

# select type of kinship estimates to use (king or pcrelate)
kin.type <- tolower(config["kinship_method"])
if (is.na(config["kinship_file"])) {
    kinMat <- divMat
} else {
    if (kin.type == "king") {
        cfg <- setNames(config["kinship_file"], "king_file")
    } else if (kin.type == "pcrelate") {
        cfg <- setNames(config["kinship_file"], "pcrelate_file")
        kinMat <- getKinship(cfg, sample.id)[[1]]
    }
    kinMat <- getKinship(cfg, sample.id)[[1]]
}
message("Using ", kin.type, " kinship estimates")

# divide into related and unrelated set
kin_thresh <- as.numeric(config["kinship_threshold"])
part <- pcairPartition(kinobj=kinMat, kin.thresh=kin_thresh,
                       divobj=divMat, div.thresh=-kin_thresh,
                       sample.include=sample.id)

rels <- part$rels
unrels <- part$unrels
save(rels, file=config["out_related_file"])
save(unrels, file=config["out_unrelated_file"])
message("Found ", length(unrels), " unrelated and ", length(rels), " related samples")

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
