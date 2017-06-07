library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(dplyr)
sessionInfo()

argp <- arg_parser("LocusZoom plots")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--segment", help="row in locus_file to plot", default=1, type="integer")
argv <- parse_args(argp)
config <- readConfig(argv$config)
segment <- argv$segment

required <- c("assoc_file",
              "locus_file")
optional <- c("flanking_region"=500,
              "out_prefix"="locuszoom")
config <- setConfigDefaults(config, required, optional)
print(config)

# read selected variant
locus <- read.table(config["locus_file"], header=TRUE, as.is=TRUE)[segment,]
variant <- locus$variantID
var.chr <- locus$chr
pop <- locus$pop
stopifnot(pop %in% c("AFR", "AMR", "ASN", "EUR"))

assocfile <- insertChromString(config["assoc_file"], var.chr)

## get association test results
assoc <- getobj(assocfile)
flank <- as.numeric(config["flanking_region"]) * 1000
var.pos <- assoc$pos[assoc$variantID == variant]
start <- var.pos - flank
end <- var.pos + flank
assoc <- assoc %>%
    filter(chr == var.chr, pos > start, pos < end) %>%
    select(variantID, chr, pos, MAF, ends_with("pval")) %>%
    mutate(MarkerName=paste0("chr", chr, ":", pos))
names(assoc)[5] <- "P-value"

## construct METAL-format file
metal <- assoc[,c("MarkerName", "P-value")]
assoc.filename <- tempfile()
write.table(metal, file=assoc.filename, row.names=FALSE, quote=FALSE, sep="\t")


ld.cmd <- paste("--pop", pop, "--source 1000G_March2012")
ld.title <- paste("1000G", pop)

prefix <- paste0(config["out_prefix"], "_var", variant, "_ld_", pop)
lz.name <- paste0("chr", var.chr, ":", var.pos)
maf <- assoc$MAF[assoc$variantID == variant]
title <- paste(lz.name, "- LD:", ld.title, "- MAF:", formatC(maf, digits=3))
#title <- paste("LD:", ld.title)

command <- paste("locuszoom",
                 "theme=publication",
                 "--cache None",
                 "--no-date",
                 "--plotonly",
                 "--gene-table gencode",
                 "--build hg19",
                 "--chr", var.chr,
                 "--metal", assoc.filename,
                 ld.cmd,
                 paste0("--refsnp \"", lz.name, "\""),
                 paste0("--flank ", config["flanking_region"], "kb"),
                 #paste("--chr", var.chr, "--start", start, "--end", end),
                 "--prefix ", prefix,
                 paste0("title=\"", title, "\""),
                 paste0("signifLine=\"", -log10(5e-8), "\" signifLineColor=\"gray\" signifLineWidth=\"2\""),
                 "ylab=\"-log10(p-value) from single variant test\"")

cat(paste(command, "\n"))
system(command)


