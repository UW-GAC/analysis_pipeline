library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(dplyr)
sessionInfo()

argp <- arg_parser("LocusZoom plots")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argp <- add_argument(argp, "--variant", help="variant ID", type="integer")
argp <- add_argument(argp, "--pop", help="population (AFR, AMR, ASN, EUR)", type="character")
argv <- parse_args(argp)
config <- readConfig(argv$config)
var.chr <- intToChr(argv$chromosome)
variant <- argv$variant
pop <- argv$pop
stopifnot(pop %in% c("AFR", "AMR", "ASN", "EUR"))

required <- c("assoc_file")
              #"gds_file")
optional <- c("flanking_region"=500,
              "out_prefix"="locuszoom")
config <- setConfigDefaults(config, required, optional)
print(config)

## gds file can have two parts split by chromosome identifier
#gdsfile <- config["gds_file"]
assocfile <- config["assoc_file"]
if (!is.na(var.chr)) {
    #gdsfile <- insertChromString(gdsfile, var.chr)
    assocfile <- insertChromString(assocfile, var.chr)
}

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

## construct EPACTS-format file
## gds <- seqOpen(gdsfile)
## seqSetFilter(gds, variant.id=assoc$variantID)
## epacts <- granges(gds) %>%
##     as.data.frame() %>%
##     cbind(ref=refChar(gds), alt=altChar(gds)) %>%
##     mutate(MARKER_ID=paste0(seqnames, ":", start, "_", ref, "/", alt)) %>%
##     left_join(assoc, by=c("seqnames"="chr", "start"="pos")) %>%
##     select(seqnames, start, end, MARKER_ID, MAF, pval, variantID) %>%
##     rename_("#CHROM"="seqnames", "BEGIN"="start", "END"="end", "PVALUE"="pval")
#seqClose(gds)

## construct METAL-format file
metal <- assoc[,c("MarkerName", "P-value")]
assoc.filename <- tempfile()
write.table(metal, file=assoc.filename, row.names=FALSE, quote=FALSE, sep="\t")


ld.cmd <- paste("--pop", pop, "--source 1000G_March2012")
ld.title <- paste("1000G", pop)

prefix <- paste0(config["out_prefix"], "_var", variant, "_ld_", pop)
#lz.name <- paste0("chr", var.chr, ":", var.pos)
#maf <- assoc$MAF[assoc$variantID == variant]
#title <- paste(lz.name, "- LD:", ld.title, "- MAF:", formatC(maf, digits=3))
title <- paste("LD:", ld.title)

command <- paste("/projects/resources/software/apps/locuszoom/bin/locuszoom",
                 "theme=publication",
                 "--cache None",
                 "--no-date",
                 "--chr", var.chr,
                 "--plotonly",
                 #"--epacts", assoc.filename,
                 "--metal", assoc.filename,
                 ld.cmd,
                 "--build hg19",
                 #paste0("--refsnp \"", lz.name, "\""),
                 #paste0("--flank ", config["flanking_region"], "kb"),
                 paste("--chr", var.chr, "--start", start, "--end", end),
                 "--prefix ", prefix,
                 paste0("title=\"", title, "\""),
                 paste0("signifLine=\"", -log10(5e-8), "\" signifLineColor=\"gray\" signifLineWidth=\"2\""))

cat(paste(command, "\n"))
system(command)


