library(TopmedPipeline)
library(SeqVarTools)
library(GENESIS)
library(gdsfmt)
sessionInfo()

args <- commandArgs(trailingOnly=TRUE)
config <- readConfig(args[1])

required <- c("gds_file",
              "pca_file",
              "variant_include_file")
optional <- c("n_pcs"=3,
              "out_prefix"="pcrelate",
              "sample_block_size"=10000,
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

pca <- getobj(config["pca_file"])

n_pcs <- min(as.integer(config["n_pcs"]), length(pca$unrels))

pcrelate(seqData,
         pcMat=pca$vectors[,1:n_pcs],
         training.set=pca$unrels,
         scan.include=sample.id, snp.include=variant.id,
         write.to.gds=TRUE, gds.prefix=config["out_prefix"],
         scan.block.size=as.integer(config["sample_block_size"]))

seqClose(seqData)
