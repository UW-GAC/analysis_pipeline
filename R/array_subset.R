library(argparser)
library(TopmedPipeline)
library(Biobase)
library(SeqVarTools)
library(GenomicRanges)
sessionInfo()

argp <- arg_parser("subset sequencing data for comparing with array")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)

required <- c("array_gds_file",
              "seq_gds_file",
              "seq_annot_file",
              "subset_gds_file",
              "study")
optional <- c("maf_threshold"=0.05,
              "missing_threshold"=0.01)
config <- setConfigDefaults(config, required, optional)
print(config)

array.gds <- seqOpen(config["array_gds_file"])

gdsfile <- config["seq_gds_file"]
outfile <- config["subset_gds_file"]
if (!is.na(chr)) {
    gdsfile <- insertChromString(gdsfile, chr)
    outfile <- insertChromString(outfile, chr, err="subset_gds_file")
}
seq.gds <- seqOpen(gdsfile)

# filter to study samples
seq.annot <- getobj(config["seq_annot_file"])
samples <- seq.annot$sample.id[seq.annot$study %in% config["study"]]
seqSetFilter(seq.gds, sample.id=samples)

# filter to overlaps with array
array.gr <- granges(array.gds)
seqSetFilter(seq.gds, variant.sel=array.gr)

# filter by MAF and missing rate
seqSetFilterCond(seq.gds,
                 maf=as.numeric(config["maf_threshold"]),
                 missing.rate=as.numeric(config["missing_threshold"]))

seqExport(seq.gds, outfile, fmt.var=character(), info.var=character())

seqClose(seq.gds)
seqClose(array.gds)
