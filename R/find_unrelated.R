library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Partition samples into related and unrelated sets")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("king_file")
optional <- c("kinship_method"="king",
              "kinship_threshold"=0.04419417, # 2^(-9/2), 3rd degree
              "out_related_file"="related.RData",
              "out_unrelated_file"="unrelated.RData",
              "pcrelate_file"=NA,
              "sample_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

## always use king estimates for ancestry divergence
king <- getobj(config["king_file"])
divMat <- king$kinship
colnames(divMat) <- rownames(divMat) <- king$sample.id

if (!is.na(config["sample_include_file"])) {
    sample.id <- as.character(getobj(config["sample_include_file"]))
    ind <- colnames(divMat) %in% sample.id
    divMat <- divMat[ind, ind]
}
message("Using ", nrow(divMat), " samples")

## select type of kinship estimates to use (king or pcrelate)
kin.type <- tolower(config["kinship_method"])
if (kin.type == "king") {
    kinMat <- divMat
} else if (kin.type == "pcrelate") {
    pcr <- openfn.gds(config["pcrelate_file"])
    kinMat <- pcrelateMakeGRM(pcr, scan.include=colnames(divMat), scaleKin=1)
    closefn.gds(pcr)
} else {
    stop("kinship method should be 'king' or 'pcrelate'")
}
message("Using ", kin.type, " kinship estimates")

# divide into related and unrelated set
kin.thresh <- as.numeric(config["kinship_threshold"])
part <- pcairPartition(kinMat=kinMat, kin.thresh=kin.thresh,
                       divMat=divMat, div.thresh=-kin.thresh)

rels <- part$rels
unrels <- part$unrels
save(rels, file=config["out_related_file"])
save(unrels, file=config["out_unrelated_file"])
message("Found ", length(unrels), " unrelated and ", length(rels), " related samples")

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
