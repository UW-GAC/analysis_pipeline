library(argparser)
library(TopmedPipeline)
library(SNPRelate)
library(GENESIS)
library(gdsfmt)
library(Biobase)
library(readr)
library(dplyr)
library(ggplot2)
library(hexbin)
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

## ## Some code to plot king ibdseg results with hexbins:
## dat <- fread('/projects/topmed/research/relatedness/analysts/caitlin/data/freeze5_round2pcairvars_allSamps_allChrs.seg')
## dat <- dat[, k0 := 1 - IBD1Seg - IBD2Seg]

## p <- ggplot(dat, aes(k0, 0.5*PropIBD)) +
##         geom_abline(intercept = 0.25, slope = -0.25) + 
##         geom_hex(aes(fill = log10(..count..)), bins = 100) +
##         scale_fill_gradientn(values = rescale(c(log10(1),log10(10),log10(1000000))), colours = c('black','steelblue','pink')) + 
##         geom_hline(yintercept = 2^(-c(3,5,7,9,11)/2), linetype = 'dashed', size = 0.3) + 
##         scale_y_continuous(breaks = c(2^(-c(2,4,6,8,10)/2),0), labels = round(c(2^(-c(2,4,6,8,10)/2),0), 3)) +
##         scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.25)) +
##         xlab('KING ibdseg k0') + ylab('KING ibdseg Kinship') + 
##         theme(panel.grid.minor = element_blank())

## select type of kinship estimates to use (king or pcrelate)
kin.type <- tolower(config["kinship_method"])
kin.thresh <- as.numeric(config["kinship_threshold"])
if (kin.type == "king_ibdseg") {
    kinship <- read_tsv(config["kinship_file"], col_types="-c-c--nnn-") %>%
        mutate(IBS0=(1 - IBD1Seg - IBD2Seg), kinship=0.5*PropIBD)
    xvar <- "IBS0"
} else if (kin.type == "king_related") {
    kinship <- read_tsv(config["kinship_file"], col_types="-cc----n--n-----") %>%
        rename(kinship=Kinship)
    xvar <- "IBS0"
} else if (kin.type == "king_kinship") {
    kinship <- read_tsv(config["kinship_file"], col_types="-cc----nn-") %>%
        rename(kinship=Kinship)
    xvar <- "IBS0"
} else if (kin.type == "king") {
    if (tools::file_ext(config["kinship_file"]) == "gds") {
        king <- gds2ibdobj(config["kinship_file"], sample.id=sample.id)
        kinship <- snpgdsIBDSelection(king, kinship.cutoff=kin.thresh)
    } else {
        king <- getobj(config["kinship_file"])
        samp.sel <- if (is.null(sample.id)) NULL else king$sample.id %in% sample.id
        kinship <- snpgdsIBDSelection(king, kinship.cutoff=kin.thresh, samp.sel=samp.sel)
    }
    xvar <- "IBS0"
    rm(king)
} else if (kin.type == "pcrelate") {
    pcr <- getobj(config["kinship_file"])
    kinship <- pcr$kinBtwn %>%
        rename(kinship=kin) %>%
        select(ID1, ID2, k0, kinship)
    xvar <- "k0"
    rm(pcr)
} else {
    stop("kinship method should be 'king' or 'pcrelate'")
}
message("Plotting ", kin.type, " kinship estimates")

p <- ggplot(kinship, aes_string(xvar, "kinship")) +
    geom_hline(yintercept=2^(-seq(3,11,2)/2), linetype="dashed", color="grey") +
    ## geom_point(alpha=0.5) +
    ## theme_bw()
    geom_hex(aes(fill = log10(..count..)), bins = 100) +
    ylab("kinship estimate")

ggsave(config["out_file_all"], plot=p, width=6, height=6)


## plot separately by study
if (!is.na(config["phenotype_file"]) & !is.na(config["study"])) {
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
        geom_hline(yintercept=2^(-seq(3,11,2)/2), linetype='dashed', color="grey") +
        ## geom_point(alpha=0.5) +
        ## theme_bw()
        facet_wrap(~study) +
        geom_hex(aes(fill = log10(..count..)), bins = 100) +
        ylab("kinship estimate")

    ggsave(config["out_file_study"], plot=p, width=7, height=7)

    ## only plot cross-study for king --ibdseg or --pcrelate
    if (kin.type %in% c("king_ibdseg", "pcrelate")) {
        # only plot cross-study relatives >= Deg2
        kinship.cross <- kinship %>%
            filter(study1 != study2) %>%
            filter(kinship > 2^(-7/2))

        # only make the plot if there are some cross-study kinship pairs
        # leave this one as geom_point instead of hexbin - color-code by study, and not many points
        if (nrow(kinship.cross) > 0){
            p <- ggplot(kinship.cross, aes_string(xvar, "kinship", color="study2")) +
                geom_hline(yintercept=2^(-seq(3,11,2)/2), linetype='dashed', color="grey") +
                geom_point() +
                facet_wrap(~study1, drop=FALSE) +
                ylab("kinship estimate") +
                theme_bw()
            ggsave(config["out_file_cross"], plot=p, width=8, height=7)
        }
    }
}

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
