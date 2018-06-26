library(argparser)
library(TopmedPipeline)
library(gdsfmt)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
sessionInfo()

argp <- arg_parser("PCA correlation plots")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("corr_file")
optional <- c("chromosomes"="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22",
              "n_pcs"=20,
              "n_perpage"=4,
              "out_prefix"="pca_corr",
              "thin"=TRUE)
config <- setConfigDefaults(config, required, optional)
print(config)

chr <- strsplit(config["chromosomes"], " ", fixed=TRUE)[[1]]
files <- sapply(chr, function(c) insertChromString(config["corr_file"], c, "corr_file"))

corr <- do.call(rbind, lapply(unname(files), function(f) {
    c <- openfn.gds(f)
    dat <- t(read.gdsn(index.gdsn(c, "correlation")))
    n_pcs <- min(as.integer(config["n_pcs"]), ncol(dat))
    dat <- dat[,1:n_pcs]
    missing <- rowSums(is.na(dat)) == n_pcs # monomorphic variants
    dat <- dat[!missing,]
    colnames(dat) <- paste0("PC", 1:n_pcs)
    dat <- data.frame(dat,
                      chr=readex.gdsn(index.gdsn(c, "chromosome"), sel=!missing),
                      pos=readex.gdsn(index.gdsn(c, "position"), sel=!missing),
                      stringsAsFactors=FALSE)
    closefn.gds(c)

    ## transform to data frame with PC as column
    dat <- dat %>%
        gather(PC, value, -chr, -pos) %>%
        filter(!is.na(value)) %>%
        mutate(value=abs(value)) %>%
        mutate(PC=factor(PC, levels=paste0("PC", 1:n_pcs)))

    ## thin points
    ## take up to 10,000 points from each of 10 evenly spaced bins
    if (as.logical(config["thin"])) {
        dat <- thinPoints(dat, "value", n=10000, nbins=10, groupBy="PC")
    }

    dat
}))

## make chromosome a factor so they are plotted in order
corr <- mutate(corr, chr=factor(chr, levels=c(1:22, "X")))
chr <- levels(corr$chr)
cmap <- setNames(rep_len(brewer.pal(8, "Dark2"), length(chr)), chr)

# plot over multiple pages
n_pcs <- length(unique(corr$PC))
n_plots <- ceiling(n_pcs/as.integer(config["n_perpage"]))
bins <- as.integer(cut(1:n_pcs, n_plots))
for (i in 1:n_plots) {
    bin <- paste0("PC", which(bins == i))
    dat <- filter(corr, PC %in% bin)

    p <- ggplot(dat, aes(chr, value, group=interaction(chr, pos), color=chr)) +
        geom_point(position=position_dodge(0.8)) +
        facet_wrap(~PC, scales="free", ncol=1) +
        scale_color_manual(values=cmap, breaks=names(cmap)) +
        ylim(0,1) +
        theme_bw() +
        theme(legend.position="none") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        xlab("chromosome") + ylab("abs(correlation)")
    ggsave(paste0(config["out_prefix"], "_" , i, ".png"), plot=p, width=10, height=15)
}

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
