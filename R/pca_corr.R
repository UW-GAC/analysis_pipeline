library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(SNPRelate)
library(gdsfmt)
sessionInfo()

argp <- arg_parser("Correlation of variants with PCs")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)

required <- c("gds_file",
              "pca_file")
optional <- c("n_pcs"=32,
              "out_file"="pca_corr.gds",
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

## if we have a chromosome indicator but only one gds file, select chromosome
if (!is.na(chr) && !bychrfile) {
    filterByChrom(gds, chr)
}

if (!is.na(varfile)) {
    filterByFile(gds, varfile)
} else {
    filterByPass(gds)
    filterBySNV(gds)
}

variant.id <- seqGetData(gds, "variant.id")
message("Using ", length(variant.id), " variants")

pca <- getobj(config["pca_file"])
n_pcs <- min(as.integer(config["n_pcs"]), length(pca$sample.id))
nt <- countThreads()
snpgdsPCACorr(pca, gdsobj=gds, snp.id=variant.id, eig.which=1:n_pcs,
              num.thread=nt, outgds=outfile)

## add chromosome and position to output
seqSetFilter(gds, variant.id=variant.id)
chromosome <- seqGetData(gds, "chromosome")
position <- seqGetData(gds, "position")
seqClose(gds)

pca.corr <- openfn.gds(outfile, readonly=FALSE)
add.gdsn(pca.corr, "chromosome", chromosome, compress="LZMA_RA")
add.gdsn(pca.corr, "position", position, compress="LZMA_RA")
closefn.gds(pca.corr)
cleanup.gds(outfile)


# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
