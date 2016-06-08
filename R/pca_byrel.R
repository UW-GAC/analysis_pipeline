library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(SNPRelate)
sessionInfo()

argp <- arg_parser("PCA on unrelated set, project relatives")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("gds_file",
              "related_file",
              "unrelated_file",
              "variant_include_file")
optional <- c("n_pcs"=20,
              "out_file"="pca.RData",
              "sample_include_file"=NA)
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

variant.id <- getobj(config["variant_include_file"])
length(variant.id)

# number of PCs to return
n_pcs <- min(as.integer(config["n_pcs"]), length(unrels))

# run PCA on unrelated set
nt <- countThreads()
pca.unrel <- snpgdsPCA(gds, sample.id=unrels, snp.id=variant.id,
                       eigen.cnt=n_pcs, num.thread=nt)

# project values for relatives
snp.load <- snpgdsPCASNPLoading(pcaobj=pca.unrel, gdsobj=gds, num.thread=nt)
samp.load <- snpgdsPCASampLoading(loadobj=snp.load, gdsobj=gds, sample.id=rels,
                                  num.thread=nt)

# combine unrelated and related PCs and order as in GDS file
eigenvect <- rbind(pca.unrel$eigenvect, samp.load$eigenvect)
rownames(eigenvect) <- c(pca.unrel$sample.id, samp.load$sample.id)
seqSetFilter(gds, sample.id=rownames(eigenvect))
sample.id <- seqGetData(gds, "sample.id")
samp.ord <- match(sample.id, rownames(eigenvect))
eigenvect <- eigenvect[samp.ord,]

# output object
pca <- list(vectors=eigenvect,
            values=pca.unrel$eigenval,
            varprop=pca.unrel$varprop,
            rels=rels,
            unrels=unrels)

save(pca, file=config["out_file"])

seqClose(gds)
