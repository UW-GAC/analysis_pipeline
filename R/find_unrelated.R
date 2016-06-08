library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Partition samples into related and unrelated sets")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("ibd_file")
optional <- c("kinship_method"="king",
              "kinship_threshold"=0.04419417, # 2^(-9/2), 3rd degree
              "n_pcs"=20,
              "out_related_file"="related.RData",
              "out_unrelated_file"="unrelated.RData",
              "pcrelate_file"=NA,
              "sample_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

ibd <- getobj(config["ibd_file"])
divMat <- ibd$kinship
colnames(divMat) <- rownames(divMat) <- ibd$sample.id

if (!is.na(config["sample_include_file"])) {
    sample.id <- as.character(getobj(config["sample_include_file"]))
    divMat <- divMat[sample.id, sample.id]
}

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

# divide into related and unrelated set
kin_thresh <- as.numeric(config["kinship_threshold"])
part <- pcairPartition(kinMat=kinMat, kin.thresh=kin_thresh,
                       divMat=divMat, div.thresh=-kin_thresh)

rels <- part$rels
unrels <- part$unrels
save(rels, file=config["out_related_file"])
save(unrels, file=config["out_unrelated_file"])

