library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(SNPRelate)
sessionInfo()

argp <- arg_parser("PCA on unrelated set, project relatives")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("gds_file",
              "related_file",
              "unrelated_file")
optional <- c("n_pcs"=32,
              "out_file"="pca.RData",
              "out_file_unrel"="pca_unrel.RData",
              "sample_include_file"=NA,
              "variant_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

gds <- seqOpen(config["gds_file"])

rels <- getobj(config["related_file"])
unrels <- getobj(config["unrelated_file"])

if (!is.na(config["sample_include_file"])) {
    sample.id <- as.character(getobj(config["sample_include_file"]))
    rels <- intersect(rels, sample.id)
    unrels <- intersect(unrels, sample.id)
}
message("Using ", length(unrels), " unrelated and ", length(rels), " related samples")

if (!is.na(config["variant_include_file"])) {
    variant.id <- getobj(config["variant_include_file"])
    message("Using ", length(variant.id), " variants")
} else {
    variant.id <- NULL
}

# number of PCs to return
n_pcs <- min(as.integer(config["n_pcs"]), length(unrels))

# run PCA on unrelated set
message("PCA on unrelated set")
nt <- countThreads()
pca.unrel <- snpgdsPCA(gds, sample.id=unrels, snp.id=variant.id,
                       algorithm="randomized", eigen.cnt=n_pcs, num.thread=nt)
save(pca.unrel, file=config["out_file_unrel"])

# project values for relatives
message("PCA projection for related set")
snp.load <- snpgdsPCASNPLoading(pca.unrel, gdsobj=gds, num.thread=nt)
samp.load <- snpgdsPCASampLoading(snp.load, gdsobj=gds, sample.id=rels, num.thread=nt)

# combine unrelated and related PCs and order as in GDS file
eigenvect <- rbind(pca.unrel$eigenvect, samp.load$eigenvect)
rownames(eigenvect) <- c(pca.unrel$sample.id, samp.load$sample.id)
seqSetFilter(gds, sample.id=rownames(eigenvect), verbose=FALSE)
sample.id <- seqGetData(gds, "sample.id")
samp.ord <- match(sample.id, rownames(eigenvect))
eigenvect <- eigenvect[samp.ord,]

# output object
pca <- list(vectors=eigenvect,
            values=pca.unrel$eigenval[1:n_pcs],
            varprop=pca.unrel$varprop[1:n_pcs],
            rels=rels,
            unrels=unrels)

save(pca, file=config["out_file"])

seqClose(gds)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
