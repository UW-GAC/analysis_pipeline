library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(GENESIS)
library(SNPRelate)
sessionInfo()

argp <- arg_parser("PC-AiR (partition relatives and run PCA using SNPRelate)")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("gds_file",
              "king_file",
              "variant_include_file")
optional <- c("kinship_threshold"=0.04419417, # 2^(-9/2), 3rd degree
              "n_pcs"=20,
              "out_file"="pca.RData",
              "sample_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

gds <- seqOpen(config["gds_file"])

king <- getobj(config["king_file"])
kinship <- king$kinship
colnames(kinship) <- rownames(kinship) <- king$sample.id

if (!is.na(config["sample_include_file"])) {
    sample.id <- as.character(getobj(config["sample_include_file"]))
    kinship <- kinship[sample.id, sample.id]
}

variant.id <- getobj(config["variant_include_file"])

# divide into related and unrelated set
kin_thresh <- as.numeric(config["kinship_threshold"])
part <- pcairPartition(kinMat=kinship, kin.thresh=kin_thresh,
                       divMat=kinship, div.thresh=-kin_thresh)

# number of PCs to return
n_pcs <- min(as.integer(config["n_pcs"]), length(part$unrels))

# run PCA on unrelated set
nt <- countThreads()
pca.unrel <- snpgdsPCA(gds, sample.id=part$unrels, snp.id=variant.id,
                       eigen.cnt=n_pcs, num.thread=nt)

# project values for relatives
snp.load <- snpgdsPCASNPLoading(pcaobj=pca.unrel, gdsobj=gds, num.thread=nt)
samp.load <- snpgdsPCASampLoading(loadobj=snp.load, gdsobj=gds, sample.id=part$rels,
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
            rels=part$rels,
            unrels=part$unrels)

save(pca, file=config["out_file"])

seqClose(gds)
