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

if (!is.null(chr)) {
    chrom <- seqGetData(gds, "chromosome")
    seqSetFilter(gds, variant.sel=(chrom == chr), verbose=FALSE)
    var.chr <- seqGetData(gds, "variant.id")
    variant.id <- intersect(variant.id, var.chr)
    seqResetFilter(gds, verbose=FALSE)
}

maf <- as.numeric(config["maf_threshold"])
r <- as.numeric(config["ld_r_threshold"])
win <- as.numeric(config["ld_win_size"]) * 1e6

snpset <- snpgdsLDpruning(gds, sample.id=sample.id, snp.id=variant.id, maf=maf, 
                          method="corr", slide.max.bp=win, ld.threshold=r,
                          num.thread=countThreads())

pruned <- unlist(snpset, use.names=FALSE)
save(pruned, file=outfile)

seqClose(gds)
