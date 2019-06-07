library(argparser)
library(TopmedPipeline)
library(Biobase)
library(SeqVarTools)
library(GenomicRanges)
sessionInfo()

argp <- arg_parser("duplicate discordance")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--segment", help="sample block", type="integer")
argp <- add_argument(argp, "--n_segments", help="number of sample blocks", type="integer")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("array_gds_file",
              "array_annot_file",
              "seq_gds_file",
              "seq_annot_file",
              "study")
optional <- c("granges_include_file"=NA,
              "sample_include_file"=NA,
              "variant_include_file"=NA,
              "out_prefix"="dup_disc")
config <- setConfigDefaults(config, required, optional)
print(config)

## sample blocks
segment <- argv$segment
n_segments <- argv$n_segments
message("sample block ", segment, " of ", n_segments)

array.gds <- seqOpen(config["array_gds_file"])
array.annot <- getobj(config["array_annot_file"])
arrayData <- SeqVarData(array.gds, sampleData=array.annot)

seq.gds <- seqOpen(config["seq_gds_file"])
seq.annot <- getobj(config["seq_annot_file"])
## allow subset GDS
seq.annot <- seq.annot[match(seqGetData(seq.gds, "sample.id"), seq.annot$sample.id),]
seqData <- SeqVarData(seq.gds, sampleData=seq.annot)

# selected variants as granges
if (!is.na(config["granges_include_file"])) {
    gr <- getobj(config["granges_include_file"])
    seqSetFilter(arrayData, variant.sel=gr)
    seqSetFilter(seqData, variant.sel=gr)
}

# for seq data, also filter on ids
if (!is.na(config["variant_include_file"])) {
    ids <- getobj(config["variant_include_file"])
    seqSetFilter(seqData, variant.id=ids, action="intersect")
}

if (is.na(config["sample_include_file"])) {
    samples <- seq.annot$sample.id[seq.annot$study %in% config["study"]]
} else {
    samples <- getobj(config["sample_include_file"])
}

# divide samples into blocks
if (n_segments > 1) {
    bins <- cut(seq_along(samples), breaks=n_segments, labels=FALSE)
    samples <- samples[bins == segment]
}

seqSetFilter(seqData, sample.id=samples)

res <- duplicateDiscordance(seqData, arrayData,
                            match.samples.on=c("submitted_subject_id", "subject.id"),
                            match.variants.on="position", discordance.type="hethom",
                            by.variant=FALSE, verbose=TRUE)

save(res, file=constructFilename(config["out_prefix"], segment=segment))
