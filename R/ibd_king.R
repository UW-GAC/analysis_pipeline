library(TopmedPipeline)
library(SeqVarTools)
library(SNPRelate)
sessionInfo()

args <- commandArgs(trailingOnly=TRUE)
config <- readConfig(args[1])

required <- c("gds_file",
              "variant_include_file")
optional <- c("out_file"="ibd_king.RData",
              "sample_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

gds <- seqOpen(config["gds_file"])

if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
} else {
    sample.id <- NULL
}

variant.id <- getobj(config["variant_include_file"])

ibd <- snpgdsIBDKING(gds, sample.id=sample.id, snp.id=variant.id,
                     num.thread=countThreads())

save(ibd, file=config["out_file"])

seqClose(gds)
