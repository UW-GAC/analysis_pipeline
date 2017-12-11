library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(GENESIS)
library(gdsfmt)
sessionInfo()

argp <- arg_parser("PC-AiR (partition relatives and run PCA)")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("gds_file",
              "king_file",
              "variant_include_file")
optional <- c("kinship_method"="king",
              "kinship_threshold"=0.04419417, # 2^(-9/2), 3rd degree
              "n_pcs"=20,
              "out_file"="pcair.RData",
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

king <- getobj(config["king_file"])
divMat <- king$kinship
colnames(divMat) <- rownames(divMat) <- king$sample.id

kin.type <- tolower(config["kinship_method"])
if (kin.type == "king") {
    kinMat <- divMat
} else if (kin.type == "pcrelate") {
    pcr <- openfn.gds(config["pcrelate_file"])
    kinMat <- pcrelateMakeGRM(pcr, scaleKin=1)
    closefn.gds(pcr)
} else {
    stop("kinship method should be 'king' or 'pcrelate'")
}
message("Using ", kin.type, " kinship estimates")

kin_thresh <- as.numeric(config["kinship_threshold"])
n_pcs <- min(as.integer(config["n_pcs"]), length(king$sample.id))

pca <- pcair(seqData, v=n_pcs,
             kinMat=kinMat, kin.thresh=kin_thresh,
             divMat=divMat, div.thresh=-kin_thresh,
             scan.include=sample.id, snp.include=variant.id)

save(pca, file=config["out_file"])

seqClose(seqData)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
