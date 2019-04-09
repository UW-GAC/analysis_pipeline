library(argparser)
library(TopmedPipeline)
library(dplyr)
library(tidyr)
library(Biobase)
library(SeqArray)
sessionInfo()

argp <- arg_parser("bed to gds")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("bed_file",
              "bim_file",
              "fam_file",
              "annot_file")
optional <- c("gds_file"="out.gds",
              "platform_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

bed.fn <- config["bed_file"]
bim.fn <- config["bim_file"]
fam.fn <- config["fam_file"]

fam <- read.table(fam.fn, header=FALSE, as.is=TRUE)
names(fam) <- c("family", "sample.id", "father", "mother", "sex", "phenotype")

# create a unique sample id
if (any(duplicated(fam$sample.id))) {
    fam <- fam %>%
        rename(id=sample.id) %>%
        group_by(id) %>%
        mutate(sample.id=paste(id, row_number(), sep="_")) %>%
        ungroup() %>%
        select(family, sample.id, father, mother, sex, phenotype)
    new.fam.fn <- paste0(fam.fn, ".unique")
    write.table(fam, new.fam.fn, quote=FALSE, row.names=FALSE, col.names=FALSE)
    fam.fn <- new.fam.fn
}

if (!is.na(config["platform_file"])) {
    platform <- read.table(config["platform_file"], header=TRUE, as.is=TRUE)

    # this will fail if subject.id contains an underscore
    fam <- fam %>%
        select(family, sample.id, sex) %>%
        separate(family, into=c("phs", "platform.id", "subject.id"), sep="_") %>%
        mutate(platform.id=as.integer(platform.id)) %>%
        left_join(platform, by=c(platform.id="ID")) %>%
        select(-platform.id)
} else {
    fam <- fam %>%
        select(family, sample.id, sex) %>%
        mutate(subject.id=sample.id)
}

annot <- AnnotatedDataFrame(fam)
save(annot, file=config["annot_file"])

seqBED2GDS(bed.fn, fam.fn, bim.fn, out.gdsfn=config["gds_file"])
