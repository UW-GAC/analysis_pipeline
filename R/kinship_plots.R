library(argparser)
library(TopmedPipeline)
library(SNPRelate)
library(GENESIS)
library(gdsfmt)
library(Biobase)
library(dplyr)
library(ggplot2)
sessionInfo()

argp <- arg_parser("Kinship plots")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("kinship_file")
optional <- c("kinship_method"="king",
              "kinship_threshold"=0.04419417, # 2^(-9/2), 3rd degree
              "out_file_all"="kinship_all.pdf",
              "out_file_cross"="kinship_cross_study.pdf",
              "out_file_study"="kinship_within_study.pdf",
              "phenotype_file"=NA,
              "sample_include_file"=NA,
              "study"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
} else {
    sample.id <- NULL
}

## select type of kinship estimates to use (king or pcrelate)
kin.type <- tolower(config["kinship_method"])
kin.thresh <- as.numeric(config["kinship_threshold"])
if (kin.type == "king") {
    ## king <- getobj(config["kinship_file"])
    ## samp.sel <- if (is.null(sample.id)) NULL else king$sample.id %in% sample.id
    ## kinship <- snpgdsIBDSelection(king, kinship.cutoff=kin.thresh, samp.sel=samp.sel)
    king <- gds2ibdobj(config["kinship_file"], sample.id=sample.id)
    kinship <- snpgdsIBDSelection(king, kinship.cutoff=kin.thresh)
    xvar <- "IBS0"
} else if (kin.type == "pcrelate") {
    pcr <- openfn.gds(config["kinship_file"])
    kinship <- pcrelateReadKinship(pcr, kin.thresh=kin.thresh, scan.include=sample.id)
    closefn.gds(pcr)
    kinship <- kinship %>%
        rename(kinship=kin) %>%
        select(ID1, ID2, k0, kinship)
    xvar <- "k0"
} else {
    stop("kinship method should be 'king' or 'pcrelate'")
}
message("Plotting ", kin.type, " kinship estimates")

p <- ggplot(kinship, aes_string(xvar, "kinship")) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
    geom_point(alpha=0.5) +
    ylab("kinship estimate") +
    theme_bw()
ggsave(config["out_file_all"], plot=p, width=6, height=6)


## plot separately by study
if (is.na(config["phenotype_file"]) | is.na(config["study"])) q("no")
study <- config["study"]
message("Plotting by study variable ", study)

annot <- getobj(config["phenotype_file"])
stopifnot(study %in% varLabels(annot))
annot <- pData(annot) %>%
    select_("sample.id", study)

kinship <- kinship %>%
    left_join(annot, by=c(ID1="sample.id")) %>%
    rename_(study1=study) %>%
    left_join(annot, by=c(ID2="sample.id")) %>%
    rename_(study2=study)

kinship.study <- kinship %>%
    filter(study1 == study2) %>%
    rename(study=study1) %>%
    select(-study2)

p <- ggplot(kinship.study, aes_string(xvar, "kinship")) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype='dashed', color="grey") +
    geom_point(alpha=0.5) +
    facet_wrap(~study) +
    ylab("kinship estimate") +
    theme_bw()
p <- ggsave(config["out_file_study"], plot=p, width=7, height=7)

kinship.cross <- kinship %>%
    filter(study1 != study2)

# only make the plot if there are some cross-study kinship pairs
if (nrow(kinship.cross) > 0){
  p <- ggplot(kinship.cross, aes_string(xvar, "kinship", color="study2")) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype='dashed', color="grey") +
    geom_point() +
    facet_wrap(~study1, drop=FALSE) +
    ylab("kinship estimate") +
    theme_bw()
  ggsave(config["out_file_cross"], plot=p, width=8, height=7)
}

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
