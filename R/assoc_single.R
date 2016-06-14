library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(Biobase)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Association test - single variant")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome number (1-24)", type="integer")
argv <- parse_args(argp)
config <- readConfig(argv$config)
chr <- argv$chromosome

required <- c("gds_file",
              "null_model_file",
              "phenotype_file")
optional <- c("mac_threshold"=30, # takes precedence
              "maf_threshold"=0.01,
              "out_file"="assoc_single.RData",
              "pass_only"=TRUE,
              "test_type"="score",
              "variant_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

## gds file can have two parts split by chromosome identifier
gdsfile <- config["gds_file"]
outfile <- config["out_file"]
if (!is.na(chr)) {
    if (chr == 23) chr <- "X"
    if (chr == 24) chr <- "Y"
    gdsfile <- insertChromString(gdsfile, chr, err="gds_file")
    outfile <- insertChromString(outfile, chr, err="out_file")
}
    
gds <- seqOpen(gdsfile)

# get null model
nullModel <- getobj(config["null_model_file"])

# get samples included in null model
sample.id <- nullModel$scanID

if (!is.na(config["variant_include_file"])) {
    variant.id <- getobj(config["variant_include_file"])
    seqSetFilter(gds, variant.id=variant.id)
    } else {
    variant.id <- seqGetData(gds, "variant.id")
}

if (as.logical(config["pass_only"])) {
    filt <- seqGetData(gds, "annotation/filter")
    variant.id <- variant.id[filt == "PASS"]
    seqSetFilter(gds, variant.id=variant.id)
}

mac.min <- as.numeric(config["mac_threshold"])
maf.min <- as.numeric(config["maf_threshold"])
if ((!is.na(mac.min) & mac.min > 1) |
    (!is.na(maf.min) & maf.min > 0)) {
    seqSetFilter(gds, sample.id=sample.id, verbose=FALSE)
    ref.freq <- seqAlleleFreq(gds)
    maf <- pmin(ref.freq, 1-ref.freq)
    if (!is.na(mac.min)) {
        maf.filt <- 2 * maf * (1-maf) * length(sample.id) >= mac.min
        message(paste("Running on", sum(maf.filt), "variants with MAC >=", mac.min))
    } else {
        maf.filt <- maf >= maf.min
        message(paste("Running on", sum(maf.filt), "variants with MAF >=", maf.min))
    }
    variant.id <- variant.id[maf.filt]
    #seqSetFilter(gds, variant.id=variant.id)
}
seqResetFilter(gds, verbose=FALSE)

# get phenotypes
annot <- getobj(config["phenotype_file"])

# createSeqVarData object
seqData <- SeqVarData(gds, sampleData=annot)

test <- switch(tolower(config["test_type"]),
               score="Score",
               wald="Wald")

assoc <- assocTestMM(seqData, nullModel, test=test,
                     snp.include=variant.id)

save(assoc, file=outfile)

seqClose(seqData)
