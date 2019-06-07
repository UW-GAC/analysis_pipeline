library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(GenomicRanges)
sessionInfo()

argp <- arg_parser("duplicate discordance")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--n", help="number of target ranges", default=10000, type="integer")
argv <- parse_args(argp)
n <- argv$n
config <- readConfig(argv$config)

required <- c("subset_gds_file",
              "granges_include_file")
optional <- c("out_prefix"="dup_disc")
config <- setConfigDefaults(config, required, optional=c())
print(config)

gds <- seqOpen(config["subset_gds_file"])
gr <- granges(gds)

fp.gr <- getobj(config["granges_include_file"])

array.fp <- subsetByOverlaps(gr, fp.gr)
message(length(array.fp), " variants overlap with fingerprints")

if (length(array.fp) < 0.9*n) {
    message("selecting additional random variants")
    keep <- sort(sample(1:length(gr), n))
    gr <- gr[keep]
    new.gr <- sort(c(fp.gr, gr))
    save(new.gr, file=paste0(config["out_prefix"], "_granges.RData"))
}

seqClose(gds)
