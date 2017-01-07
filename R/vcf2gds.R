library(argparser)
library(TopmedPipeline)
library(SeqArray)
library(tools)
sessionInfo()

argp <- arg_parser("Convert VCF to GDS")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argv <- parse_args(argp)
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)

required <- c("vcf_file", "gds_file")
optional <- c(format="GT")
config <- setConfigDefaults(config, required, optional)
print(config)

## vcf file can have two parts split by chromosome identifier
vcffile <- config["vcf_file"]
gdsfile <- config["gds_file"]
if (!is.na(chr)) {
    vcffile <- insertChromString(vcffile, chr, "vcf_file")
    gdsfile <- insertChromString(gdsfile, chr, "gds_file")
}

## pick format fields to import
fmt.import <- strsplit(config["format"], " ", fixed=TRUE)[[1]]

## write to the scratch disk of each node
gdsfile.tmp <- tempfile()
message("gds temporarily located at ", gdsfile.tmp)

## is this a bcf file?
isBCF <- file_ext(vcffile) == "bcf"
if (isBCF) {
    ## use bcftools to read text
    vcffile <- pipe(paste("bcftools view", vcffile), "rt")
}

seqVCF2GDS(vcffile, gdsfile.tmp, fmt.import=fmt.import, storage.option="LZMA_RA",
           parallel=countThreads())

if (isBCF) close(vcffile)

## copy it
file.copy(gdsfile.tmp, gdsfile)
## remove the tmp file
file.remove(gdsfile.tmp)
