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
optional <- c("mac_threshold"=5, # takes precedence
              "maf_threshold"=0.001,
              "out_file"="assoc_single.RData",
              "pass_only"=TRUE,
              "test_type"="score",
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
    
gds <- seqOpen(gdsfile)

# get null model
nullModel <- getobj(config["null_model_file"])

# get samples included in null model
sample.id <- nullModel$scanID

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

variant.id <- seqGetData(gds, "variant.id")
seqResetFilter(gds, verbose=FALSE)


# get phenotypes
annot <- getobj(config["phenotype_file"])

# createSeqVarData object
seqData <- SeqVarData(gds, sampleData=annot)

test <- switch(tolower(config["test_type"]),
               score="Score",
               wald="Wald")

assoc <- assocTestMM(seqData, nullModel, test=test, snp.include=variant.id)

## make output consistent with aggregate tests
names(assoc)[names(assoc) == "snpID"] <- "variantID"
names(assoc) <- sub(".Stat", ".stat", names(assoc), fixed=TRUE)
seqSetFilter(seqData, variant.id=assoc$variantID, verbose=FALSE)
assoc$pos <- seqGetData(seqData, "position")
cols <- setdiff(names(assoc), c("variantID", "chr", "pos"))
assoc <- assoc[,c("variantID", "chr", "pos", cols)]

save(assoc, file=outfile)

seqClose(seqData)
