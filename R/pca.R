library(TopmedPipeline)
library(SeqVarTools)
library(SNPRelate)
sessionInfo()

args <- commandArgs(trailingOnly=TRUE)
config <- readConfig(args[1])

required <- c("gds_file",
              "variant_include_file")
optional <- c("n_pcs"=20,
              "out_file"="pca.RData",
              "sample_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

gds <- seqOpen(config["gds_file"])

if (!is.na(config["sample_include_file"])) {
    sample.id <- as.character(getobj(config["sample_include_file"]))
} else {
    sample.id <- NULL
}

variant.id <- getobj(config["variant_include_file"])

# run PCA on unrelated set
n_pcs <- min(as.integer(config["n_pcs"]), length(sample.id))
nt <- countThreads()
pca <- snpgdsPCA(gds, sample.id=sample.id, snp.id=variant.id,
                 eigen.cnt=n_pcs, num.thread=nt)

save(pca, file=config["out_file"])

seqClose(gds)
