library(argparser)
library(TopmedPipeline)
library(Biobase)
library(dplyr)
library(ggplot2)
library(GGally)
sessionInfo()

argp <- arg_parser("PCA plots")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("pca_file")
optional <- c("n_pairs"=6,
              "out_file_scree"="pca_scree.pdf",
              "out_file_pc12"="pca_pc12.pdf",
              "out_file_parcoord"="pca_parcoord.pdf",
              "out_file_pairs"="pca_pairs.png",
              "phenotype_file"=NA,
              "group"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)


## get PCs
pca <- getobj(config["pca_file"])
pcs <- as.data.frame(pca$vectors[pca$unrels,])
n <- ncol(pcs)
names(pcs) <- paste0("PC", 1:n)
pcs$sample.id <- row.names(pcs)

## scree plot
dat <- data.frame(pc=1:n, varprop=pca$varprop)
p <- ggplot(dat, aes(x=factor(pc), y=100*varprop)) +
  geom_point() + theme_bw() +
  xlab("PC") + ylab("Percent of variance accounted for")
ggsave(config["out_file_scree"], plot=p, width=6, height=6)

## color by group
if (!is.na(config["phenotype_file"]) & !is.na(config["group"])) {
    group <- config["group"]
    annot <- getobj(config["phenotype_file"])
    stopifnot(group %in% varLabels(annot))
    annot <- pData(annot) %>%
        select_("sample.id", group)
    pcs <- left_join(pcs, annot, by="sample.id")
} else {
    ## make up dummy group
    group <- "group"
    pcs$group <- "NA"
}

p <- ggplot(pcs, aes_string("PC1", "PC2", color=group)) + geom_point(alpha=0.5) +
    guides(colour=guide_legend(override.aes=list(alpha=1)))
ggsave(config["out_file_pc12"], plot=p, width=7, height=6)


npr <- min(as.integer(config["n_pairs"]), n)
p <- ggpairs(pcs, mapping=aes_string(color=group), columns=1:npr,
             lower=list(continuous=wrap("points", alpha=0.5)),
             diag=list(continuous="densityDiag"),
             upper=list(continuous="blank"))
png(config["out_file_pairs"], width=8, height=8, units="in", res=150)
print(p)
dev.off()


pc2 <- pcs
names(pc2)[1:ncol(pc2)] <- sub("PC", "", names(pc2)[1:ncol(pc2)])

p <- ggparcoord(pc2, columns=1:n, groupColumn=group, alphaLines=0.5, scale="uniminmax") +
    guides(colour=guide_legend(override.aes=list(alpha=1, size=2))) +
    xlab("PC") + ylab("")
ggsave(config["out_file_parcoord"], plot=p, width=10, height=5)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
