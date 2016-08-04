library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(Biobase)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Association test - sliding window")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome number (1-24)", type="integer")
argv <- parse_args(argp)
config <- readConfig(argv$config)
chr <- argv$chromosome

# add parameters for:
# user-specified weights

required <- c("gds_file",
              "null_model_file",
              "phenotype_file")
optional <- c("alt_freq_range"="0 1",
              "out_file"="assoc_window.RData",
              "pass_only"=TRUE,
              "pval_skat"="kuonen",
              "rho"="0",
              "test"="burden",
              "test_type"="score",
              "variant_include_file"=NA,
              "weight_beta"="0.5 0.5",
              "window_size"=50,
              "window_step"=20)
config <- setConfigDefaults(config, required, optional)
print(config)

## gds file can have two parts split by chromosome identifier
gdsfile <- config["gds_file"]
outfile <- config["out_file"]
varfile <- config["variant_include_file"]
if (!is.na(chr)) {
    if (chr == 23) chr <- "X"
    if (chr == 24) chr <- "Y"
    bychrfile <- grepl(" ", gdsfile) # do we have one file per chromosome?
    gdsfile <- insertChromString(gdsfile, chr)
    outfile <- insertChromString(outfile, chr, err="out_file")
    varfile <- insertChromString(varfile, chr)
}
    
gds <- seqOpen(gdsfile)

# get null model
nullModel <- getobj(config["null_model_file"])

# get samples included in null model
sample.id <- nullModel$scanID

if (!is.na(varfile)) {
    variant.id <- getobj(varfile)
    seqSetFilter(gds, variant.id=variant.id)
    } else {
    variant.id <- seqGetData(gds, "variant.id")
}

## if we have a chromosome indicator but only one gds file, select chromosome
if (!is.na(chr) & !bychrfile) {
    chrom <- seqGetData(gds, "chromosome")
    variant.id <- variant.id[chrom == chr]
    seqSetFilter(gds, variant.id=variant.id)
}

if (as.logical(config["pass_only"])) {
    filt <- seqGetData(gds, "annotation/filter")
    variant.id <- variant.id[filt == "PASS"]
    #seqSetFilter(gds, variant.id=variant.id)
}

message("Using ", length(variant.id), " variants")
seqResetFilter(gds, verbose=FALSE)


# get phenotypes
annot <- getobj(config["phenotype_file"])

# createSeqVarData object
seqData <- SeqVarData(gds, sampleData=annot)

test <- switch(tolower(config["test"]),
               burden="Burden",
               skat="SKAT")

test.type <- switch(tolower(config["test_type"]),
                    firth="Firth",
                    score="Score",
                    wald="Wald")

af.range <- as.numeric(strsplit(config["alt_freq_range"], " ", fixed=TRUE)[[1]])
weights <- as.numeric(strsplit(config["weight_beta"], " ", fixed=TRUE)[[1]])
rho <- as.numeric(strsplit(config["rho"], " ", fixed=TRUE)[[1]])
pval <- tolower(config["pval_skat"])
size <- as.numeric(config["window_size"])
step <- as.numeric(config["window_step"])

assoc <- assocTestSeqWindow(seqData, nullModel,
                            test=test,
                            burden.test=test.type,
                            AF.range=af.range,
                            weight.beta=weights,
                            rho=rho,
                            pval.method=pval,
                            variant.include=variant.id,
                            window.size=size,
                            window.shift=step)

save(assoc, file=outfile)

seqClose(seqData)
