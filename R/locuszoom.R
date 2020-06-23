library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(dplyr)
sessionInfo()

argp <- arg_parser("LocusZoom plots")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--segment", help="row in locus_file to plot", default=1, type="integer")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
segment <- argv$segment

required <- c("assoc_file",
              "locus_file")
optional <- c("flanking_region"=500,
              "gds_file"=NA,
              "genome_build"="hg38",
              "ld_sample_include"=NA,
              "locus_type"="variant",
              "out_prefix"="locuszoom",
              "signif_line"=5e-8,
              "track_file"=NA,
              "track_file_type"="window",
              "track_label"="",
              "track_threshold"=5e-8)
config <- setConfigDefaults(config, required, optional)
print(config)

stopifnot(config["locus_type"] %in% c("variant", "region"))

# read selected locus
locus <- read.table(config["locus_file"], header=TRUE, as.is=TRUE)[segment,]
stopifnot(all(c("chr", "pop") %in% names(locus)))
print(locus)

# population for LD
pop <- toupper(locus$pop)
stopifnot(pop %in% c("TOPMED", "AFR", "AMR", "ASN", "EUR", "EAS", "SAS"))

## get association test results
var.chr <- locus$chr
assocfile <- insertChromString(config["assoc_file"], var.chr)
assoc <- getobj(assocfile)

if (config["locus_type"] == "variant") {
    stopifnot("variant.id" %in% names(locus))
    variant <- locus$variant.id
    flank <- as.numeric(config["flanking_region"]) * 1000
    var.pos <- assoc$pos[assoc$variant.id == variant]
    start <- var.pos - flank
    if (start < 1) start <- 1
    end <- var.pos + flank
    
    lz.name <- paste0("chr", var.chr, ":", var.pos)
    ld.region <- paste0("--refsnp \"", lz.name, "\"", " --flank ", config["flanking_region"], "kb")
    prefix <- paste0(config["out_prefix"], "_var", variant, "_ld_", pop)
    freq <- assoc$freq[assoc$variant.id == variant]
    maf <- min(freq, 1-freq)
    mac <- assoc$MAC[assoc$variant.id == variant]
    title <- paste(lz.name, "- MAF:", formatC(maf, digits=3), "- MAC:", mac)
    
} else if (config["locus_type"] == "region") {
    stopifnot(all(c("start", "end") %in% names(locus)))
    start <- locus$start
    end <- locus$end

    ld.region <- paste("--chr", var.chr, "--start", start, "--end", end)
    prefix <- paste0(config["out_prefix"], "_ld_", pop)
    title <- ""
}

## construct METAL-format file
assoc <- assoc %>%
    filter(chr == var.chr, pos > start, pos < end) %>%
    select(variant.id, chr, pos, ends_with("pval"))
names(assoc)[4] <- "pval"

##### NOTE: i added the following two lines of code #####
## remove duplicate chr:pos rows, removing the variant with the less significant (higher) pvalue
assoc <- assoc[order(assoc$pval,decreasing=FALSE),]
assoc <- assoc[!duplicated(assoc$pos),]
assoc <- assoc[order(assoc$pos),]
#####

assoc.filename <- tempfile()
writeMETAL(assoc, file=assoc.filename)

# LD
if (pop != "TOPMED") {
    ld.cmd <- paste("--pop", pop, "--source 1000G_Nov2014")
    ld.title <- paste("LD: 1000G", pop)
} else {
    if (!is.na(config["ld_sample_include"])) {
        sample.id <- getobj(config["ld_sample_include"])
    } else {
        sample.id <- NULL
    }
    if (config["locus_type"] == "variant") {
        ref.var <- variant
    } else {
        ref.var <- filter(assoc, pval == min(pval))$variant.id
        if (length(ref.var) > 1) {
            message("Multiple variants with minimum pval; selecting the first one as reference")
            ref.var <- ref.var[1]
        }
    }
    gdsfile <- insertChromString(config["gds_file"], var.chr)
    ld <- calculateLD(gdsfile, variant.id=assoc$variant.id, ref.var=ref.var, sample.id=sample.id)
    ld.filename <- tempfile()
    writeLD(assoc, ld, ref.var, file=ld.filename)

    ld.cmd <- paste("--ld", ld.filename)
    ld.title <- "LD: TOPMed"
}
title <- if (title == "") ld.title else paste(ld.title, title, sep=" - ")

## construct BED track file
if (!is.na(config["track_file"])) {
    trackfile <- insertChromString(config["track_file"], var.chr)
    track <- getAssoc(trackfile, config["track_file_type"]) %>%
        filter(pval < as.numeric(config["track_threshold"]))
    track.filename <- tempfile()
    writeBED(track, file=track.filename, track.label=config["track_label"])
    track.cmd <- paste("--bed-tracks", track.filename)
} else {
    track.cmd <- ""
}

signif <- as.numeric(config["signif_line"])

command <- paste("locuszoom",
                 "theme=publication",
                 "--cache None",
                 "--no-date",
                 "--plotonly",
                 "--gene-table gencode",
                 "--build", config["genome_build"],
                 "--chr", var.chr,
                 "--metal", assoc.filename,
                 track.cmd,
                 ld.cmd,
                 ld.region,
                 "--prefix ", prefix,
                 paste0("title=\"", title, "\""),
                 paste0("signifLine=\"", -log10(signif), "\" signifLineColor=\"gray\" signifLineWidth=\"2\""),
                 "ylab=\"-log10(p-value) from single variant test\"")

cat(paste(command, "\n"))
system(command)

unlink(assoc.filename)
if (exists("track.filename")) unlink(track.filename)
if (exists("ld.filename")) unlink(ld.filename)
