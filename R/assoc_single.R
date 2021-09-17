library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(Biobase)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Association test - single variant")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argp <- add_argument(argp, "--segment", help="segment number", type="integer")
argp <- add_argument(argp, "--num_cores", help="number of cores", type="integer", default=1)
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)
segment <- argv$segment

# parallelization
if (argv$num_cores > 1) {
    BPPARAM <- BiocParallel::MulticoreParam(workers=argv$num_cores)
} else {
    BPPARAM <- BiocParallel::SerialParam()
}

required <- c("gds_file",
              "null_model_file",
              "phenotype_file")
optional <- c("genome_build"="hg38",
              "mac_threshold"=5, # takes precedence
              "maf_threshold"=0.001,
              "out_prefix"="assoc_single",
              "pass_only"=TRUE,
              "segment_file"=NA,
              "test_type"="score",
              "variant_include_file"=NA,
              "variant_block_size"=1024,
              "genotype_coding"="additive")
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(basename(argv$config), ".assoc_single.params"))

## gds file can have two parts split by chromosome identifier
gdsfile <- config["gds_file"]
varfile <- config["variant_include_file"]
if (!is.na(chr)) {
    bychrfile <- grepl(" ", gdsfile) # do we have one file per chromosome?
    gdsfile <- insertChromString(gdsfile, chr)
    varfile <- insertChromString(varfile, chr)
}

gds <- seqOpen(gdsfile)

# get phenotypes
annot <- getobj(config["phenotype_file"])

# createSeqVarData object
annot <- matchAnnotGds(gds, annot)
seqData <- SeqVarData(gds, sampleData=annot)

# get null model
nullModel <- getobj(config["null_model_file"])

# get samples included in null model
nullModel <- GENESIS:::.updateNullModelFormat(nullModel)
sample.id <- nullModel$fit$sample.id

# check if null model is for fast.score.SE
if(isNullModelFastScore(nullModel)){
  fast.score.SE <- TRUE
  message("Using fast approximation to score standard error")
}else{
  fast.score.SE <- FALSE
}

if (!is.na(segment)) {
    filterBySegment(seqData, segment, config["segment_file"])
}

if (!is.na(varfile)) {
    filterByFile(seqData, varfile)
}

## if we have a chromosome indicator but only one gds file, select chromosome
if (!is.na(chr) && !bychrfile) {
    filterByChrom(seqData, chr)
}

if (as.logical(config["pass_only"])) {
    filterByPass(seqData)
}

## MAC/MAF filtering
mac.min <- as.numeric(config["mac_threshold"])
maf.min <- as.numeric(config["maf_threshold"])
build <- config["genome_build"]
if (!is.na(mac.min)) {
    filterByMAC(seqData, sample.id, mac.min=mac.min, build=build)
} else {
    filterByMAF(seqData, sample.id, maf.min=maf.min, build=build)
}

checkSelectedVariants(seqData)

# create iterator
block.size <- as.integer(config["variant_block_size"])
iterator <- SeqVarBlockIterator(seqData, variantBlock=block.size)

test <- switch(tolower(config["test_type"]),
               score="Score",
               score.spa="Score.SPA",
               binomirare="BinomiRare")

geno.coding <- config["genotype_coding"]

assoc <- assocTestSingle(iterator, nullModel,
                         test=test,
                         fast.score.SE=fast.score.SE,
                         genome.build=build,
                         BPPARAM=BPPARAM,
                         geno.coding=geno.coding)

save(assoc, file=constructFilename(config["out_prefix"], chr, segment))

seqClose(seqData)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
