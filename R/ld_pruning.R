library(TopmedPipeline)
library(SeqVarTools)
library(SNPRelate)
sessionInfo()

args <- commandArgs(trailingOnly=TRUE)
config <- readConfig(args[1])

required <- c("gds_file")
optional <- c("ld_r_threshold"=0.32,
              "ld_win_size"=10,
              "maf_threshold"=0.01,
              "out_file"="pruned_variants.RData",
              "sample_include_file"=NA,
              "variant_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

gds <- seqOpen(config["gds_file"])

if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
} else {
    sample.id <- NULL
}

if (!is.na(config["variant_include_file"])) {
    variant.id <- getobj(config["variant_include_file"])
} else {
    filt <- seqGetData(gds, "annotation/filter") == "PASS"
    snv <- isSNV(gds, biallelic=TRUE)
    variant.id <- seqGetData(gds, "variant.id")[filt & snv]
}

maf <- as.numeric(config["maf_threshold"])
r <- as.numeric(config["ld_r_threshold"])
win <- as.numeric(config["ld_win_size"]) * 1e6

snpset <- snpgdsLDpruning(gds, sample.id=sample.id, snp.id=variant.id, maf=maf, 
                          method="corr", slide.max.bp=win, ld.threshold=r,
                          num.thread=countThreads())

pruned <- unlist(snpset, use.names=FALSE)
save(pruned, file=config["out_file"])

seqClose(gds)
