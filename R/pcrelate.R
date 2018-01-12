library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(GENESIS)
library(gdsfmt)
sessionInfo()

argp <- arg_parser("PC-Relate")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

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
    message("Using ", length(sample.id), " samples")
} else {
    sample.id <- NULL
    message("Using all samples")
}

variant.id <- getobj(config["variant_include_file"])
message("Using ", length(variant.id), " variants")

pca <- getobj(config["pca_file"])

n_pcs <- min(as.integer(config["n_pcs"]), length(pca$unrels))

pcrelate(seqData,
         pcMat=pca$vectors[,1:n_pcs],
         training.set=pca$unrels,
         scan.include=sample.id, snp.include=variant.id,
         write.to.gds=TRUE, gds.prefix=config["out_prefix"],
         scan.block.size=as.integer(config["sample_block_size"]))

seqClose(seqData)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
