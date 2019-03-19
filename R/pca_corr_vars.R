library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
sessionInfo()

argp <- arg_parser("Select variants for correlation with PCs")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)

required <- c("gds_file",
              "segment_file")
optional <- c("n_corr_vars"=10e6,
              "out_file"="pca_corr_variants.RData",
              "variant_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

## gds file can have two parts split by chromosome identifier
gdsfile <- config["gds_file"]
outfile <- config["out_file"]
varfile <- config["variant_include_file"]
if (!is.na(chr)) {
    message("Running on chromosome ", chr)
    bychrfile <- grepl(" ", gdsfile) # do we have one file per chromosome?
    gdsfile <- insertChromString(gdsfile, chr)
    outfile <- insertChromString(outfile, chr, err="out_file")
    varfile <- insertChromString(varfile, chr)
}

gds <- seqOpen(gdsfile)

filterByPass(gds)
filterBySNV(gds)

vars <- seqGetData(gds, "variant.id")

seqClose(gds)

## get number of variants with same proportion of this chromosome in genome
segments <- read.table(config["segment_file"], header=TRUE, stringsAsFactors=FALSE)
prop <- sum(segments$chromosome == chr)/nrow(segments)
nvars <- round(as.integer(config["n_corr_vars"]) * prop)

if (length(vars) > nvars) {
    vars <- sort(sample(vars, nvars))
}
message("selected ", length(vars), " variants")

## add pruned variants
if (!is.na(varfile)) {
    pruned <- getobj(varfile)
    vars <- sort(unique(c(vars, pruned)))
}
message(length(vars), " total variants")


save(vars, file=outfile)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
