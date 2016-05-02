library(TopmedPipeline)
library(SeqVarTools)
library(GENESIS)
sessionInfo()

args <- commandArgs(trailingOnly=TRUE)
config <- readConfig(args[1])

required <- c("gds_file",
              "ibd_file",
              "variant_include_file")
optional <- c("kinship_threshold"=0.04419417, # 2^(-9/2), 3rd degree
              "n_pcs"=20,
              "out_file"="pcair.RData",
              "sample_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

gds <- seqOpen(config["gds_file"])
seqData <- SeqVarData(gds)

if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
} else {
    sample.id <- NULL
}

variant.id <- getobj(config["variant_include_file"])

ibd <- getobj(config["ibd_file"])
kinship <- ibd$kinship
colnames(kinship) <- rownames(kinship) <- ibd$sample.id

kin_thresh <- as.numeric(config["kinship_threshold"])
n_pcs <- min(as.integer(config["n_pcs"]), length(ibd$sample.id))

pca <- pcair(seqData, v=n_pcs, kin.thresh=kin_thresh,
             kinMat=kinship, divMat=kinship,
             scan.include=sample.id, snp.include=variant.id)

save(pca, file=config["out_file"])
    
seqClose(seqData)
