library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(GENESIS)
sessionInfo()

argp <- arg_parser("PC-Relate")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("gds_file",
              "pca_file")
optional <- c("n_pcs"=3,
              "out_prefix"="pcrelate_beta",
              "sample_include_file"=NA,
              "variant_block_size"=1024,
              "variant_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

gds <- seqOpen(config["gds_file"])
seqData <- SeqVarData(gds)

if (!is.na(config["variant_include_file"])) {
    filterByFile(seqData, config["variant_include_file"])
}

pca <- getobj(config["pca_file"])
n_pcs <- min(as.integer(config["n_pcs"]), length(pca$unrels))
pcs <- as.matrix(pca$vectors[,1:n_pcs])
sample.include <- samplesGdsOrder(seqData, pca$unrels)

if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
    sample.include <- intersect(sample.include, sample.id)
}


# create iterator
block.size <- as.integer(config["variant_block_size"])
iterator <- SeqVarBlockIterator(seqData, variantBlock=block.size)

beta <- calcISAFBeta(iterator,
                     pcs=pcs,
                     sample.include=sample.include)

save(beta, file=paste0(config["out_prefix"], ".RData"))

seqClose(seqData)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
