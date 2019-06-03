library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(SNPRelate)
sessionInfo()

argp <- arg_parser("IBD with KING-robust")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("gds_file")
optional <- c("out_file"="ibd_king.gds",
              #"out_file"="ibd_king.RData",
              "sample_include_file"=NA,
              "variant_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

gds <- seqOpen(config["gds_file"])

if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
    message("Using ", length(sample.id), " samples")
} else {
    sample.id <- NULL
    message("Using all samples")
}

if (!is.na(config["variant_include_file"])) {
    variant.id <- getobj(config["variant_include_file"])
    message("Using ", length(variant.id), " variants")
} else {
    variant.id <- NULL
    message("Using all variants")
}

ibd <- snpgdsIBDKING(gds, sample.id=sample.id, snp.id=variant.id,
                     num.thread=countThreads())

#save(ibd, file=config["out_file"])
list2gds(ibd, config["out_file"])

seqClose(gds)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
