library(argparser)
library(TopmedPipeline)
library(Biobase)
library(SeqVarTools)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Calculate variant score and SE")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-22)", type="character")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)

required <- c("gds_file",
			  "null_model_file",
			  "phenotype_file")
optional <- c("genome_build"="hg38",
			  "n_variants"=10,
			  "mac_threshold"=20,
			  "pass_only"=TRUE,
			  "variant_include_file"=NA,
			  "out_prefix"="variant_score_table")
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(basename(argv$config), ".variant_score_calc.params"))

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

checkSelectedVariants(seqData)


# user specified variants (if provided)
if(!is.na(varfile)){
	variant.id <- seqGetData(seqData, 'variant.id')
}else{
	variant.id <- NULL
}

# compute se.ratio 
# if is.na(varfile): for a random subset of variants
# if !is.na(varfile): for specified variants
tab <- calcScore(seqData, 
				 nullModel,
				 variant.id = variant.id,
				 nvar = config["n_variants"],
				 min.mac = config["mac_threshold"],
				 genome.build = config["genome_build"])

save(tab, file = constructFilename(config["out_prefix"], chr, NA))

seqClose(seqData)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
