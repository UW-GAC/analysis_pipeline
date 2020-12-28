library(argparser)
library(TopmedPipeline)
library(Biobase)
library(SeqVarTools)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Update null model for fast score approximation association tests")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("gds_file",
			  "null_model_file",
			  "phenotype_file")
optional <- c("genome_build"="hg38",
			  "n_variants"=10,
			  "mac_threshold"=20,
			  "pass_only"=TRUE,
			  "out_prefix"="null_model")
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(basename(argv$config), ".null_model_fast_scoreSE.params"))

# get null model
nullModel <- getobj(config["null_model_file"])

# compute se.ratio for a random subset of variants from each chr
tab <- NULL
# loop through chromosomes
for(chr in 1:22){
	# gds file
	gdsfile <- config["gds_file"]
	gdsfile <- insertChromString(gdsfile, chr)
	gds <- seqOpen(gdsfile)

	# get phenotypes
	annot <- getobj(config["phenotype_file"])

	# createSeqVarData object
	annot <- matchAnnotGds(gds, annot)
	seqData <- SeqVarData(gds, sampleData=annot)

	if (as.logical(config["pass_only"])) {
	    filterByPass(seqData)
	}

	checkSelectedVariants(seqData)

	chr.tab <- calcScore(seqData, 
						 nullModel, 
						 nvar = config["n_variants"],
						 min.mac = config["mac_threshold"],
						 genome.build = config["genome_build"])
	tab <- rbind(tab, chr.tab)

	seqClose(seqData)
}

# update the null model
nullModel <- nullModelFastScore(nullModel, tab)
outfile <- sprintf("%s_fast_scoreSE.RData", config["out_prefix"])
save(nullModel, file = outfile)


# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
