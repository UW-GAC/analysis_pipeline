library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
sessionInfo()

argp <- arg_parser("Parse table of variants of regions to a list")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)

required <- c("variant_group_file")
optional <- c("aggregate_type"="allele",
              "group_id"="group_id",
              "out_file"="aggregate_list.RData")
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(basename(argv$config), ".aggregate_list.params"))

## file can have two parts split by chromosome identifier
outfile <- config["out_file"]
varfile <- config["variant_group_file"]
if (!is.na(chr)) {
    outfile <- insertChromString(outfile, chr, err="out_file")
    varfile <- insertChromString(varfile, chr)
}

groups <- getobj(varfile)

## rename columns if necessary
names(groups)[names(groups) %in% config["group_id"]] <- "group_id"
names(groups)[names(groups) %in% c("chromosome", "CHROM")] <- "chr"
names(groups)[names(groups) %in% c("position", "POS")] <- "pos"
names(groups)[names(groups) %in% "REF"] <- "ref"
names(groups)[names(groups) %in% "ALT"] <- "alt"

## subset groups by chromosome
groups <- groups[groups$chr == chr,]

if (!is.character(groups$chr)) groups$chr <- as.character(groups$chr)

if (config["aggregate_type"] == "allele") {
    aggVarList <- aggregateGRangesList(groups)
} else if (config["aggregate_type"] == "position") {
    aggVarList <- aggregateGRanges(groups)
} else {
    stop("aggregrate_type must be 'allele' or 'position'")
}

save(aggVarList, file=outfile)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
