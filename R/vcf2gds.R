library(TopmedPipeline)
library(SeqArray)
sessionInfo()

args <- commandArgs(trailingOnly=TRUE)
config <- readConfig(args[1])

required <- c("vcf_file", "gds_file")
optional <- c()
config <- setConfigDefaults(config, required, optional)
print(config)

## is this an array job by chromosome?
chr <- if (length(args) > 1) args[2] else NULL

## vcf file can have two parts split by chromosome identifier
vcffile <- config["vcf_file"]
gdsfile <- config["gds_file"]
if (!is.null(chr)) {
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
