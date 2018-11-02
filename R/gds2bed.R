library(argparser)
library(TopmedPipeline)
library(SeqArray)
library(SNPRelate)
sessionInfo()

argp <- arg_parser("Convert GDS to BED")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("gds_file", "bed_file")
optional <- c("variant_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

gds <- seqOpen(config["gds_file"])

if (!is.na(config["variant_include_file"])) {
    filterByFile(seqData, config["variant_include_file"])
}

snpfile <- tempfile()
seqGDS2SNP(gds, snpfile)
seqClose(gds)

gds <- snpgdsOpen(snpfile)
snpgdsGDS2BED(gds, config["bed_file"])
snpgdsClose(gds)

unlink(snpfile)
