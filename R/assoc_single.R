library(TopmedPipeline)
library(SeqVarTools)
library(Biobase)
library(GENESIS)
sessionInfo()

args <- commandArgs(trailingOnly=TRUE)
config <- readConfig(args[1])

required <- c("gds_file",
              "null_model_file")
optional <- c("maf_threshold"=0.01,
              "out_file"="assoc_single.RData",
              "pass_only"=TRUE,
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
    gdsfile <- insertChromString(gdsfile, chr, "gds_file")
    outfile <- insertChromString(outfile, chr, "out_file")
}
    
gds <- seqOpen(config["gds_file"])

# get null model
nullModel <- getobj(config["null_model_file"])

# get samples included in null model
sample.id <- nullModel$scanID

if (!is.na(config["variant_include_file"])) {
    variant.id <- getobj(config["variant_include_file"])
    } else {
    variant.id <- seqGetData(gds, "variant.id")
}

if (as.logical(config["pass_only"])) {
    filt <- seqGetData(gds, "annotation/filter")
    seqSetFilter(gds, variant.sel=(filt == "PASS"), verbose=FALSE)
    var.filt <- seqGetData(gds, "variant.id")
    variant.id <- intersect(variant.id, var.filt)
    seqResetFilter(gds, verbose=FALSE)
}

maf.min <- as.numeric(config["maf_threshold"])
if (maf.min > 0) {
    seqSetFilter(gds, variant.id=variant.id, sample.id=sample.id, verbose=FALSE)
    ref.freq <- seqAlleleFreq(gds)
    maf <- pmin(ref.freq, 1-ref.freq)
    variant.id <- variant.id[maf >= maf.min]
    seqResetFilter(gds, verbose=FALSE)
}


# get phenotypes
annot <- getobj(config["phenotype_file"])

# createSeqVarData object
seqData <- SeqVarData(gds, sampleData=annot)

test <- if (nullModel$family$family == "gaussian") "Wald" else "Score"

assoc <- assocTestMM(seqData, nullModel, test=test, snp.include=variant.id)

save(assoc, file=config["out_file"])

seqClose(seqData)
