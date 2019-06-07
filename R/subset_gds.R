library(argparser)
library(TopmedPipeline)
library(SeqArray)
sessionInfo()

argp <- arg_parser("Subset GDS file")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)

required <- c("gds_file",
              "subset_gds_file")
optional <- c("sample_include_file"=NA,
              "variant_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

## gds file can have two parts split by chromosome identifier
gdsfile <- config["gds_file"]
outfile <- config["subset_gds_file"]
varfile <- config["variant_include_file"]
if (!is.na(chr)) {
    gdsfile <- insertChromString(gdsfile, chr)
    outfile <- insertChromString(outfile, chr, err="subset_gds_file")
    varfile <- insertChromString(varfile, chr)
}

gds <- seqOpen(gdsfile)

if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
} else {
    sample.id <- NULL
}

if (!is.na(varfile)) {
    variant.id <- getobj(varfile)
} else {
    variant.id <- NULL
}

seqSetFilter(gds, sample.id=sample.id, variant.id=variant.id)
seqExport(gds, outfile, fmt.var=character(), info.var=character())

seqClose(gds)
