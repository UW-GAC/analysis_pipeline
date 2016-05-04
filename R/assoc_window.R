library(TopmedPipeline)
library(SeqVarTools)
library(Biobase)
library(GENESIS)
sessionInfo()

args <- commandArgs(trailingOnly=TRUE)
config <- readConfig(args[1])

# add parameters for:
# window size and step
# type of test for burden
# user-specified weights

required <- c("gds_file",
              "null_model_file",
              "phenotype_file")
optional <- c("alt_freq_range"="0 1",
              "out_file"="assoc_window.RData",
              "pass_only"=TRUE,
              "rho"="0",
              "test"="burden",
              "variant_include_file"=NA,
              "weights"="0.5 0.5")
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

# get null model
nullModel <- getobj(config["null_model_file"])

# get samples included in null model
sample.id <- nullModel$scanID

if (!is.na(config["variant_include_file"])) {
    variant.id <- getobj(config["variant_include_file"])
    } else {
    variant.id <- seqGetData(gds, "variant.id")
}

if (!is.null(chr)) {
    chrom <- seqGetData(gds, "chromosome")
    seqSetFilter(gds, variant.sel=(chrom == chr), verbose=FALSE)
    var.chr <- seqGetData(gds, "variant.id")
    variant.id <- intersect(variant.id, var.chr)
    seqResetFilter(gds, verbose=FALSE)
}

if (as.logical(config["pass_only"])) {
    filt <- seqGetData(gds, "annotation/filter")
    seqSetFilter(gds, variant.sel=(filt == "PASS"), verbose=FALSE)
    var.filt <- seqGetData(gds, "variant.id")
    variant.id <- intersect(variant.id, var.filt)
    seqResetFilter(gds, verbose=FALSE)
}


# get phenotypes
annot <- getobj(config["phenotype_file"])

# createSeqVarData object
seqData <- SeqVarData(gds, sampleData=annot)

test <- switch(tolower(config["test"]),
               burden="Burden",
               skat="SKAT")

af.range <- as.numeric(strsplit(config["alt_freq_range"], " ", fixed=TRUE)[[1]])
weights <- as.numeric(strsplit(config["weights"], " ", fixed=TRUE)[[1]])
rho <- as.numeric(strsplit(config["rho"], " ", fixed=TRUE)[[1]])

assoc <- assocTestSeqWindow(seqData, nullModel, test=test,
                            AF.range=af.range,
                            weight.beta=weights, rho=rho,
                            variant.include=variant.id)

save(assoc, file=outfile)

seqClose(seqData)
