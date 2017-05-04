library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(Biobase)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Association test - sliding window")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argp <- add_argument(argp, "--segment", help="segment number", type="integer")
argv <- parse_args(argp)
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)
segment <- argv$segment

# add parameters for:
# user-specified weights

required <- c("gds_file",
              "null_model_file",
              "phenotype_file")
optional <- c("alt_freq_range"="0 1",
              "out_prefix"="assoc_window",
              "pass_only"=TRUE,
              "pval_skat"="kuonen",
              "rho"="0",
              "segment_file"=NA,
              "test"="burden",
              "test_type"="score",
              "variant_include_file"=NA,
              "weight_beta"="0.5 0.5",
              "window_size"=50,
              "window_step"=20)
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(argv$config, "assoc_window.params"))

## gds file can have two parts split by chromosome identifier
gdsfile <- config["gds_file"]
varfile <- config["variant_include_file"]
if (!is.na(chr)) {
    bychrfile <- grepl(" ", gdsfile) # do we have one file per chromosome?
    gdsfile <- insertChromString(gdsfile, chr)
    varfile <- insertChromString(varfile, chr)
}
    
gds <- seqOpen(gdsfile)

# get null model
nullModel <- getobj(config["null_model_file"])

# get samples included in null model
sample.id <- nullModel$scanID

size <- as.numeric(config["window_size"])
step <- as.numeric(config["window_step"])

if (!is.na(segment)) {
    ## pad each segment by window size to be sure we get all possible windows
    filterBySegment(gds, segment, config["segment_file"], pad.right=size*1000)
}

if (!is.na(varfile)) {
    filterByFile(gds, varfile)
}

## if we have a chromosome indicator but only one gds file, select chromosome
if (!is.na(chr) && !bychrfile) {
    filterByChrom(gds, chr)
}

if (as.logical(config["pass_only"])) {
    filterByPass(gds)
}

checkSelectedVariants(gds)
variant.id <- seqGetData(gds, "variant.id")
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

save(assoc, file=constructFilename(config["out_prefix"], chr, segment))

seqClose(seqData)
