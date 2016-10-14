library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(Biobase)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Association test - single variant, unrelated")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome number (1-24)", type="integer")
argv <- parse_args(argp)
config <- readConfig(argv$config)
chr <- argv$chromosome

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
              "out_file"="assoc_single_unrel.RData",
              "pass_only"=TRUE,
              "sample_include_file"=NA,
              "test_type"="linear",
              "variant_include_file"=NA)
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

if (!is.na(varfile)) {
    variant.id <- getobj(varfile)
    seqSetFilter(gds, variant.id=variant.id)
}

## if we have a chromosome indicator but only one gds file, select chromosome
if (!is.na(chr) && !bychrfile) {
    gds <- filterByChrom(gds, chr)
}

if (as.logical(config["pass_only"])) {
    gds <- filterByPass(gds)
}

mac.min <- as.numeric(config["mac_threshold"])
maf.min <- as.numeric(config["maf_threshold"])
gds <- filterByMAF(gds, sample.id, mac.min, maf.min)

seqSetFilter(gds, sample.id=sample.id)

# createSeqVarData object
seqData <- SeqVarData(gds, sampleData=annot)

assoc <- regression(seqData, outcome=outcome, covar=covars, model.type=test)

## make output consistent with mixed model
names(assoc)[names(assoc) == "variant.id"] <- "variantID"
names(assoc) <- sub(".Stat", ".stat", names(assoc), fixed=TRUE)
names(assoc) <- sub(".Pval", ".pval", names(assoc), fixed=TRUE)
seqSetFilter(seqData, variant.id=assoc$variantID, verbose=FALSE)
assoc$chr <- seqGetData(seqData, "chromosome")
assoc$pos <- seqGetData(seqData, "position")
assoc$MAF <- pmin(assoc$freq, 1 - assoc$freq)
assoc$minor.allele <- ifelse(assoc$freq > 0.5, "ref", "alt")
init.cols <- c("variantID", "chr", "pos", "n", "MAF", "minor.allele")
cols <- setdiff(names(assoc), c(init.cols, "freq"))
assoc <- assoc[,c(init.cols, cols)]

save(assoc, file=outfile)

seqClose(seqData)
