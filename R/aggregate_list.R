library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
sessionInfo()

argp <- arg_parser("Parse table of variants of regions to a list")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argv <- parse_args(argp)
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)

required <- c("gds_file",
              "variant_group_file")
optional <- c("aggregate_type"="allele",
              "out_file"="aggregate_list.RData")
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(argv$config, ".aggregate_list.params"))

## gds file can have two parts split by chromosome identifier
gdsfile <- config["gds_file"]
outfile <- config["out_file"]
varfile <- config["variant_group_file"]
if (!is.na(chr)) {
    bychrfile <- grepl(" ", gdsfile) # do we have one file per chromosome?
    gdsfile <- insertChromString(gdsfile, chr)
    outfile <- insertChromString(outfile, chr, err="out_file")
    varfile <- insertChromString(varfile, chr)
}
    
gds <- seqOpen(gdsfile)

groups <- getobj(varfile)

## subset groups by chromosome
groups <- groups[groups$chromosome == chr,]

## chromosome must be character to match with gds
if (!is.character(groups$chromosome)) groups$chromosome <- as.character(groups$chromosome)

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
