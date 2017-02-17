library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(SNPRelate)
sessionInfo()

argp <- arg_parser("GRM")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("gds_file",
              "variant_include_file")
optional <- c("maf_threshold"=0.01,
              "method"="gcta",
              "out_file"="grm.RData",
              "sample_include_file"=NA)
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

variant.id <- getobj(config["variant_include_file"])
message("Using ", length(variant.id), " variants")

maf.min <- as.numeric(config["maf_threshold"])

method <- switch(tolower(config["method"]),
                 gcta="GCTA",
                 eigmix="EIGMIX",
                 beta="W&Z15")

grm <- snpgdsGRM(gds, sample.id=sample.id, snp.id=variant.id,
                 maf=maf.min, method=method,
                 num.thread=countThreads())

save(grm, file=config["out_file"])

seqClose(gds)
