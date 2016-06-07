library(argparser)
library(TopmedPipeline)
library(SeqArray)
sessionInfo()

argp <- arg_parser("Convert VCF to GDS")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome number (1-24)", type="integer")
argv <- parse_args(argp)
config <- readConfig(argv$config)
chr <- argv$chromosome

required <- c("vcf_file", "gds_file")
optional <- c()
config <- setConfigDefaults(config, required, optional)
print(config)

## vcf file can have two parts split by chromosome identifier
vcffile <- config["vcf_file"]
gdsfile <- config["gds_file"]
if (!is.na(chr)) {
    if (chr == 23) chr <- "X"
    if (chr == 24) chr <- "Y"
    vcffile <- insertChromString(vcffile, chr, "vcf_file")
    gdsfile <- insertChromString(gdsfile, chr, "gds_file")
}

## write to the scratch disk of each node
gdsfile.tmp <- tempfile()
message("gds temporarily located at ", gdsfile.tmp)

## import genotype only
seqVCF2GDS(vcffile, gdsfile.tmp, fmt.import="GT", parallel=countThreads())

## copy it
file.copy(gdsfile.tmp, gdsfile)
## remove the tmp file
file.remove(gdsfile.tmp)
