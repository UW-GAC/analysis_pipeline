library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(SNPRelate)
sessionInfo()

argp <- arg_parser("GRM")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)

required <- c("gds_file")
optional <- c("exclude_pca_corr"=TRUE,
              "genome_build"="hg38",
              "maf_threshold"=0.001,
              "missing_threshold"=0.01,
              "method"="gcta",
              "out_file"="grm.RData",
              "sample_include_file"=NA,
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

if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
    message("Using ", length(sample.id), " samples")
} else {
    sample.id <- NULL
    message("Using all samples")
}

if (!is.na(varfile)) {
    filterByFile(gds, varfile)
}

## if we have a chromosome indicator but only one gds file, select chromosome
if (!is.na(chr) && !bychrfile) {
    filterByChrom(gds, chr)
}

filterByPass(gds)
filterBySNV(gds)
if (as.logical(config["exclude_pca_corr"])) {
    filterByPCAcorr(gds, build=config["genome_build"])
}

variant.id <- seqGetData(gds, "variant.id")
message("Using ", length(variant.id), " variants")

maf.min <- as.numeric(config["maf_threshold"])
miss <- as.numeric(config["missing_threshold"])

method <- switch(tolower(config["method"]),
                 gcta="GCTA",
                 eigmix="EIGMIX",
                 indivbeta="IndivBeta")

## write to the scratch disk of each node
outfile.tmp <- tempfile()
message("gds temporarily located at ", outfile.tmp)

snpgdsGRM(gds, sample.id=sample.id, snp.id=variant.id,
          maf=maf.min, missing.rate=miss,
          method=method, out.fn=outfile.tmp,
          autosome.only=FALSE,
          num.thread=countThreads())

seqClose(gds)

## copy it
file.copy(outfile.tmp, outfile)
## remove the tmp file
file.remove(outfile.tmp)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
