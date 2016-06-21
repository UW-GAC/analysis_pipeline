library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
sessionInfo()

argp <- arg_parser("Parse table of variants of regions to a list")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome number (1-24)", type="integer")
argv <- parse_args(argp)
config <- readConfig(argv$config)
chr <- argv$chromosome

required <- c("gds_file",
              "variant_group_file")
optional <- c("aggregate_type"="allele",
              "out_file"="aggregate_list.RData")
config <- setConfigDefaults(config, required, optional)
print(config)

## gds file can have two parts split by chromosome identifier
gdsfile <- config["gds_file"]
outfile <- config["out_file"]
varfile <- config["variant_group_file"]
if (!is.na(chr)) {
    if (chr == 23) chr <- "X"
    if (chr == 24) chr <- "Y"
    gdsfile <- insertChromString(gdsfile, chr, err="gds_file")
    outfile <- insertChromString(outfile, chr, err="out_file")
    varfile <- insertChromString(varfile, chr, err="variant_group_file")
}
    
gds <- seqOpen(gdsfile)

groups <- getobj(varfile)

if (config["aggregate_type"] == "allele") {
    message("Sorting ", nrow(groups), " variant alleles")
    aggVarList <- aggregateListByAllele(gds, groups)
} else if (config["aggregate_type"] == "position") {
    message("Finding variants in ", nrow(groups), " groups")
    aggVarList <- aggregateListByPosition(gds, groups)
} else {
    stop("aggregrate_type must be 'allele' or 'position'")
}

save(aggVarList, file=outfile)
message("Saved ", sum(sapply(aggVarList, nrow)), " variant alleles in ", length(aggVarList), " groups")

seqClose(gds)
