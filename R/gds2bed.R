library(argparser)
library(TopmedPipeline)
library(SeqArray)
library(SNPRelate)
sessionInfo()

argp <- arg_parser("Convert GDS to BED")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("gds_file", "bed_file")
optional <- c("sample_include_file"=NA,
              "variant_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

gds <- seqOpen(config["gds_file"])

if (!is.na(config["variant_include_file"])) {
    filterByFile(gds, config["variant_include_file"])
}

if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
    seqSetFilter(gds, sample.id=sample.id)
}

snpfile <- tempfile()
seqGDS2SNP(gds, snpfile)
seqClose(gds)

gds <- snpgdsOpen(snpfile)
snpgdsGDS2BED(gds, config["bed_file"])
snpgdsClose(gds)

unlink(snpfile)
