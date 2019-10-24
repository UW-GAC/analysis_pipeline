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
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)
segment <- argv$segment

required <- c("gds_file",
              "null_model_file",
              "phenotype_file")
optional <- c("alt_freq_max"=1,
              "genome_build"="hg38",
              "out_prefix"="assoc_window",
              "pass_only"=TRUE,
              "rho"="0",
              "segment_file"=NA,
              "test"="burden",
              "variant_include_file"=NA,
              "variant_weight_file"=NA,
              "weight_user"="weight",
              "weight_beta"="1 1",
              "window_size"=50,
              "window_step"=20)
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(basename(argv$config), ".assoc_window.params"))

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

# get weights
if (!is.na(config["variant_weight_file"])) {
    dat <- getobj(config["variant_weight_file"])
    weight.user <- config["weight_user"]
    stopifnot(weight.user %in% names(dat))
    seqData <- addVariantData(seqData, dat)
    rm(dat)
} else {
    weight.user <- NULL
}

# get null model
nullModel <- getobj(config["null_model_file"])

# get samples included in null model
sample.id <- nullModel$sample.id

size <- as.numeric(config["window_size"])*1000
step <- as.numeric(config["window_step"])*1000

if (!is.na(segment)) {
    ## pad each segment by window size to be sure we get all possible windows
    filterBySegment(seqData, segment, config["segment_file"], pad.right=size)
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

af.max <- as.numeric(config["alt_freq_max"])
build <- config["genome_build"]
filterByRare(seqData, sample.id, af.max=af.max, build=build)

checkSelectedVariants(seqData)


iterator <- SeqVarWindowIterator(seqData, windowSize=size, windowShift=step)

test <- switch(tolower(config["test"]),
               burden="Burden",
               skat="SKAT",
               smmat="SMMAT",
               fastskat="fastSKAT",
               fastsmmat="fastSMMAT",
               skato="SKATO")

weight.beta <- as.numeric(strsplit(config["weight_beta"], " ", fixed=TRUE)[[1]])
rho <- as.numeric(strsplit(config["rho"], " ", fixed=TRUE)[[1]])

assoc <- assocTestAggregate(iterator, nullModel,
                            AF.max=af.max,
                            weight.beta=weight.beta,
                            weight.user=weight.user,
                            test=test,
                            rho=rho,
                            genome.build=build)

assoc <- addMAC(assoc, "window")

save(assoc, file=constructFilename(config["out_prefix"], chr, segment))

seqClose(seqData)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
