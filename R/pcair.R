library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(GENESIS)
library(gdsfmt)
sessionInfo()

argp <- arg_parser("PC-AiR (partition relatives and run PCA)")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("gds_file",
              "king_file",
              "variant_include_file")
optional <- c("kinship_method"="king",
              "kinship_threshold"=0.04419417, # 2^(-9/2), 3rd degree
              "n_pcs"=20,
              "out_file"="pcair.RData",
              "out_related_file"="related.RData",
              "out_unrelated_file"="unrelated.RData",
              "pcrelate_file"=NA,
              "sample_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

gds <- seqOpen(config["gds_file"])
seqData <- SeqVarData(gds)

if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
    message("Using ", length(sample.id), " samples")
} else {
    sample.id <- NULL
    message("Using all samples")
}

variant.id <- getobj(config["variant_include_file"])
message("Using ", length(variant.id), " variants")

# getKinship returns a list
divMat <- getKinship(config["king_file"], sample.id)[[1]]

kin.type <- tolower(config["kinship_method"])
if (kin.type == "king") {
    kinMat <- divMat
} else if (kin.type == "pcrelate") {
    kinMat <- getKinship(config["pcrelate"], sample.id)[[1]]
} else {
    stop("kinship method should be 'king' or 'pcrelate'")
}
message("Using ", kin.type, " kinship estimates")

kin_thresh <- as.numeric(config["kinship_threshold"])
n_pcs <- min(as.integer(config["n_pcs"]), nrow(divMat))

pca <- pcair(seqData,
             kinobj=kinMat, kin.thresh=kin_thresh,
             divobj=divMat, div.thresh=-kin_thresh,
             sample.include=sample.id, snp.include=variant.id,
             eigen.cnt=n_pcs)

save(pca, file=config["out_file"])

rels <- pca$rels
unrels <- pca$unrels
save(rels, file=config["out_related_file"])
save(unrels, file=config["out_unrelated_file"])

seqClose(seqData)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
