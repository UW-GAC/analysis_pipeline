library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(Biobase)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Association test - aggregate")
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
              "phenotype_file",
              "aggregate_variant_file")
optional <- c("alt_freq_range"="0 1",
              "out_prefix"="assoc_aggregate",
              "pval_skat"="kuonen",
              "rho"="0",
              "segment_file"=NA,
              "test"="burden",
              "test_type"="score",
              "variant_include_file"=NA,
              "weight_beta"="1 1")
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(basename(argv$config), ".assoc_aggregate.params"))

## gds file can have two parts split by chromosome identifier
gdsfile <- config["gds_file"]
aggfile <- config["aggregate_variant_file"]
varfile <- config["variant_include_file"]
if (!is.na(chr)) {
    gdsfile <- insertChromString(gdsfile, chr)
    aggfile <- insertChromString(aggfile, chr, err="aggregate_variant_file")
    varfile <- insertChromString(varfile, chr)
}

gds <- seqOpen(gdsfile)

# get null model
nullModel <- getobj(config["null_model_file"])

# get samples included in null model
sample.id <- nullModel$scanID

# get aggregate list
aggVarList <- getobj(aggfile)

# keep units that start in the requested segment
if (!is.na(segment)) {
    aggVarList <- subsetBySegment(aggVarList, segment, config["segment_file"])
}
if (length(aggVarList) == 0) {
    message("No aggregate units selected. Exiting gracefully.")
    q(save="no", status=0)
}

# subset to included variants
if (!is.na(varfile)) {
    variant.id <- getobj(varfile)
    aggVarList <- lapply(aggVarList, function(x) x[x$variant.id %in% variant.id,])
}

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

assoc <- assocTestSeq(seqData, nullModel, aggVarList,
                      test=test,
                      burden.test=test.type,
                      AF.range=af.range,
                      weight.beta=weights,
                      rho=rho,
                      pval.method=pval)

save(assoc, file=constructFilename(config["out_prefix"], chr, segment))

seqClose(seqData)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
