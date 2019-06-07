library(argparser)
library(TopmedPipeline)
library(Biobase)
library(dplyr)
library(SeqVarTools)
library(GenomicRanges)
sessionInfo()

argp <- arg_parser("look for a sample match among all array samples")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--segment", help="sample block", type="integer")
argp <- add_argument(argp, "--subset_variants", help="subset to 10,000 variants?", flag=TRUE)
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("array_gds_file",
              "array_annot_file",
              "seq_gds_file",
              "seq_annot_file",
              "sample_file")
optional <- c("granges_include_file"=NA,
              "variant_include_file"=NA,
              "out_prefix"="dup_disc")
config <- setConfigDefaults(config, required, optional)
print(config)


array.gds <- seqOpen(config["array_gds_file"])
array.annot <- getobj(config["array_annot_file"])
arrayData <- SeqVarData(array.gds, sampleData=array.annot)

seq.gds <- seqOpen(config["seq_gds_file"])
seq.annot <- getobj(config["seq_annot_file"])
## allow subset GDS
seq.annot <- seq.annot[match(seqGetData(seq.gds, "sample.id"), seq.annot$sample.id),]
seqData <- SeqVarData(seq.gds, sampleData=seq.annot)

# selected variants as granges
if (!is.na(config["granges_include_file"])) {
    gr <- getobj(config["granges_include_file"])
    seqSetFilter(arrayData, variant.sel=gr)
    seqSetFilter(seqData, variant.sel=gr)
}

# for seq data, also filter on ids
if (!is.na(config["variant_include_file"])) {
    ids <- getobj(config["variant_include_file"])
    seqSetFilter(seqData, variant.id=ids, action="intersect")
}

# select sample to match
segment <- argv$segment
samples <- getobj(config["sample_file"])
sample <- samples[segment]
seqSetFilter(seqData, sample.id=sample)



gr1 <- SeqVarTools:::.getGRanges(seqData)
gr2 <- SeqVarTools:::.getGRanges(arrayData)
overlappingVariants <- SeqVarTools:::.matchVariants(gr1, gr2, match.on="position")

# should we subset the number of variants to avoid running out of memory?
if (argv$subset_variants) {
    ind <- sample(1:nrow(overlappingVariants), 10000)
    overlappingVariants <- overlappingVariants[ind,]
}

# prepare genotype data frame
seqSetFilter(seqData, variant.id=overlappingVariants$variant.id.1)
dos1 <- refDosage(seqData)
seqSetFilter(arrayData, variant.id=overlappingVariants$variant.id.2)
dos2 <- refDosage(arrayData)

# order genotypes appropriately
dos1 <- dos1[, as.character(overlappingVariants$variant.id.1)]
dos2 <- dos2[, as.character(overlappingVariants$variant.id.2)]

# recode if the alt/ref alleles are switched
dos2[,overlappingVariants$recode] <- 2 - dos2[,overlappingVariants$recode]

class1 <- SeqVarTools:::.getGenotypeClass(dos1)

# remove missing genotypes
n.variants <- rep(NA,nrow(dos2))
n.concordant <- rep(NA,nrow(dos2))
for (i in 1:nrow(dos2)) {
    sel <- !is.na(dos2[i,])
    class2 <- SeqVarTools:::.getGenotypeClass(dos2[i,])

    n.variants[i] <- sum(sel)
    n.concordant[i] <- sum(SeqVarTools:::.getMatchesHetHom(class1, class2)[sel])
}


res <- pData(sampleData(arrayData)) %>%
    rename_(array.id="sample.id") %>%
    mutate(sample.id=sample) %>%
    cbind(n.variants, n.concordant)

save(res, file=paste0(config["out_prefix"], "_", sample, ".RData"))
