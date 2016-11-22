library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(Biobase)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Association test - single variant, unrelated")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome number (1-24)", type="integer")
argp <- add_argument(argp, "--segment", help="segment number", type="integer")
argv <- parse_args(argp)
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)
segment <- argv$segment

required <- c("gds_file",
              "outcome",
              "pca_file",
              "phenotype_file")
optional <- c("binary"=FALSE,
              "covars"=NA,
              "inverse_normal"=FALSE,
              "mac_threshold"=5, # takes precedence
              "maf_threshold"=0.001,
              "n_pcs"=3,
              "out_prefix"="assoc_single_unrel",
              "pass_only"=TRUE,
              "sample_include_file"=NA,
              "segment_file"=NA,
              "test_type"="linear",
              "variant_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

## gds file can have two parts split by chromosome identifier
gdsfile <- config["gds_file"]
varfile <- config["variant_include_file"]
if (!is.na(chr)) {
    bychrfile <- grepl(" ", gdsfile) # do we have one file per chromosome?
    gdsfile <- insertChromString(gdsfile, chr)
    varfile <- insertChromString(varfile, chr)
}

# get phenotypes
phen <- getPhenotypes(config)
annot <- phen[["annot"]]
outcome <- phen[["outcome"]]
covars <- phen[["covars"]]
sample.id <- phen[["sample.id"]]
message("Model: ", outcome, " ~ ", paste(covars, collapse=" + "))

test <- tolower(config["test_type"])

if (as.logical(config["binary"])) {
    stopifnot(all(annot[[outcome]] %in% c(0,1)))
    stopifnot(test %in% c("logistic", "firth"))
    family <- binomial
} else {
    stopifnot(test == "linear")
    family <- gaussian
}

if (as.logical(config["inverse_normal"])) {
    nullmod <- fitNullReg(annot, outcome=outcome, covars=covars,
                          scan.include=sample.id, family=family)
    annot <- addInvNorm(annot, nullmod, outcome, covars)
    outcome <- "resid.norm"
    covars <- NULL
}
    
gds <- seqOpen(gdsfile)

if (!is.na(segment)) {
    filterBySegment(gds, segment, config["segment_file"])
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

mac.min <- as.numeric(config["mac_threshold"])
maf.min <- as.numeric(config["maf_threshold"])
filterByMAF(gds, sample.id, mac.min, maf.min)

checkSelectedVariants(gds)

seqSetFilter(gds, sample.id=sample.id)

# createSeqVarData object
seqData <- SeqVarData(gds, sampleData=annot)

assoc <- regression(seqData, outcome=outcome, covar=covars, model.type=test)

## make output consistent with mixed model
assoc <- formatAssocSingle(seqData, assoc)

save(assoc, file=constructFilename(config["out_prefix"], chr, segment))

seqClose(seqData)
