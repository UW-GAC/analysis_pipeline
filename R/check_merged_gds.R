library(argparser)
library(TopmedPipeline)
library(SeqArray)
library(tools)
sessionInfo()

argp <- arg_parser("Check merged GDS")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)

required <- c("gds_file", "merged_gds_file")
config <- setConfigDefaults(config, required, optional=NULL)
print(config)

## gds file has two parts split by chromosome identifier
gdsfile <- config["gds_file"]
if (!is.na(chr)) {
    gdsfile <- insertChromString(gdsfile, chr, "gds_file")
}

## the chromosome-specific GDS file
gds <- seqOpen(gdsfile)
message("File: ", basename(gdsfile), "\n")
hash_chr <- seqDigest(gds, "genotype")
message("MD5: ", hash_chr, "\n")
seqClose(gds)

## the merged GDS file
gds_merged_fn <- config["merged_gds_file"]
gds_merge <- seqOpen(gds_merged_fn)
message("File: ", basename(gds_merged_fn), "\n")
seqSetFilterChrom(gds_merge, chr)
hash_merge <- seqDigest(gds_merge, "genotype")
message("MD5: ", hash_merge, "\n")
seqClose(gds_merge)


# check
stopifnot(hash_merge == hash_chr)

message("MD5 checking [OK]\n")
