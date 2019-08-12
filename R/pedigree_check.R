library(argparser)
library(TopmedPipeline)
library(Biobase)
library(SNPRelate)
library(readr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
sessionInfo()


argp <- arg_parser("Pedigree check")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("phenotype_file",
              "exp_rel_file",
              "kinship_file")
optional <- c("subjectID"="submitted_subject_id",
              "kinship_method"="king_ibdseg",
              "out_file"="kinship_obsrel.RData",
              "out_plot"="kinship.pdf",
              "sample_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

annot <- getobj(config["phenotype_file"])
annot <- pData(annot) %>%
    select(sample.id, Individ=!!config["subjectID"])

if (!is.na(config["sample_include_file"])) {
    ids <- getobj(config["sample_include_file"])
    annot <- filter(annot, sample.id %in% ids)
} else {
    ids <- annot$sample.id
}

## read ibd file
kin.type <- tolower(config["kinship_method"])
if (kin.type == "pcrelate") {
    pcr <- getobj(config["kinship_file"])
    ibd <- pcr$kinBtwn %>%
        filter(ID1 %in% ids & ID2 %in% ids) %>%
        rename(kinship=kin) %>%
        select(ID1, ID2, k0, kinship) %>%
        mutate(obs.rel=GWASTools::ibdAssignRelatednessKing(k0, kinship, cut.ibs0.err=0.1))
    xvar <- "k0"
    rm(pcr)
} else if (kin.type == "king" ) {
    king <- getobj(config["kinship_file"])
    ibd <- snpgdsIBDSelection(king, samp.sel=which(king$sample.id %in% ids)) %>%
        mutate(obs.rel=GWASTools::ibdAssignRelatednessKing(IBS0, kinship))
    xvar <- "IBS0"
    rm(king)
} else if (kin.type == "king_ibdseg") {
    ibd <- read_tsv(config["kinship_file"], col_types="-c-c--nnnc") %>%
        filter(ID1 %in% ids & ID2 %in% ids) %>% 
        mutate(IBS0=(1 - IBD1Seg - IBD2Seg), kinship=0.5*PropIBD) %>%
        mutate(obs.rel=recode(InfType, "Dup/MZ" = "Dup", "2nd" = "Deg2", "3rd" = "Deg3", "4th" = "U", "UN" = "U")) %>%
        select(ID1, ID2, IBS0, kinship, obs.rel)
    xvar <- "IBS0"
}
nrow(ibd)

## expected relatives
rel <- getobj(config["exp_rel_file"])
rel <- rel %>%
    mutate(pair=GWASTools::pasteSorted(Individ1, Individ2)) %>%
    select(one_of("pair", "Individ1", "Individ2", "family", "relation", "exp.rel", "MZtwinID"))

ibd <- select(ibd, "ID1", "ID2", !!xvar, "kinship", "obs.rel") %>%
    left_join(annot, by=c(ID1="sample.id")) %>%
    rename(Individ1=Individ) %>%
    left_join(annot, by=c(ID2="sample.id")) %>%
    rename(Individ2=Individ) %>%
    mutate(pair=GWASTools::pasteSorted(Individ1, Individ2)) %>%
    left_join(select(rel, -starts_with("Individ")), by="pair")

unobs <- rel %>%
    inner_join(annot, by=c(Individ1="Individ")) %>%
    rename(ID1=sample.id) %>%
    inner_join(annot, by=c(Individ2="Individ")) %>%
    rename(ID2=sample.id) %>%
    filter(!(pair %in% ibd$pair)) %>%
    mutate(obs.rel="U")

ibd <- bind_rows(ibd, unobs) %>%
    select(-pair) %>%
    mutate(exp.rel=ifelse(is.na(exp.rel), "U", exp.rel),
           exp.rel=ifelse(Individ1 == Individ2, "Dup", exp.rel)) %>%
    filter(!(exp.rel == "U" & obs.rel == "U"))

nrow(ibd)
save(ibd, file=config["out_file"])

table(ibd$exp.rel, ibd$obs.rel)

## plot unexpected relatives
ibd <- mutate(ibd, unexp=ifelse(exp.rel == obs.rel, "expected", "unexpected"))

rels <- c("Dup", "PO", "FS", "Deg1", "Deg2", "Deg3", "Q", "U")
cols <- c(brewer.pal(length(rels)-1, "Dark2")[c(1, 2, 3, 6, 5, 4, 7)], "black")
cmap <- setNames(cols, rels)

theme_set(theme_bw() + theme(legend.position=c(1, 1), legend.justification=c(1,1), legend.background = element_rect(colour = "black")))

p <- ggplot(ibd, aes_string(xvar, "kinship", color="exp.rel")) + facet_wrap(~unexp) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype='dashed', color="grey") +
    geom_point(alpha=0.7) +
    scale_color_manual(values=cmap, breaks=names(cmap), na.value="grey30") +
    guides(colour=guide_legend(override.aes=list(alpha=1))) +
    ylab("kinship estimate")

ggsave(config["out_plot"], plot=p, width=12, height=6)

