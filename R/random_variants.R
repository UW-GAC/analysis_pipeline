library(TopmedPipeline)
library(SeqVarTools)
library(SNPRelate)
sessionInfo()

args <- commandArgs(trailingOnly=TRUE)
config <- readConfig(args[1])

required <- c("gds_file")
optional <- c("maf_threshold"=0.05,
              "n_variants"=200000,
              "out_file"="random_variants.RData")
config <- setConfigDefaults(config, required, optional)
print(config)

## is this an array job by chromosome?
chr <- if (length(args) > 1) args[2] else NULL

## gds file can have two parts split by chromosome identifier
gdsfile <- config["gds_file"]
outfile <- config["out_file"]
if (!is.null(chr)) {
    if (chr == 23) chr <- "X"
    if (chr == 24) chr <- "Y"
    gdsfile <- insertChromString(gdsfile, chr)
    outfile <- insertChromString(outfile, chr, err="out_file")
}
    
gds <- seqOpen(gdsfile)

filt <- seqGetData(gds, "annotation/filter") == "PASS"
snv <- isSNV(gds, biallelic=TRUE)
variant.sel <- filt & snv

nvar <- as.integer(config["n_variants"])
if (!is.null(chr)) {
    chrom <- seqGetData(gds, "chromosome")
    variant.sel <- variant.sel & (chrom == chr)
    nvar <- ceiling(nvar * (sum(chrom == chr)/length(chrom)))
}

seqSetFilter(gds, variant.sel=variant.sel)
variant.id <- seqGetData(gds, "variant.id")

maf.min <- as.numeric(config["maf_threshold"])
if (maf.min > 0) {
    ref.freq <- seqAlleleFreq(gds)
    maf <- pmin(ref.freq, 1-ref.freq)
    variant.id <- variant.id[maf > maf.min]
}

if (length(variant.id) > nvar) {
    variant.id <- sample(variant.id, nvar)
}

save(variant.id, file=outfile)

seqClose(gds)
