library(argparser)
library(TopmedPipeline)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
sessionInfo()

argp <- arg_parser("PCA correlation plots")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("corr_file")
optional <- c("chromosomes"="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22",
              "n_pcs"=20,
              "n_perpage"=4,
              "out_prefix"="pca_corr")
config <- setConfigDefaults(config, required, optional)
print(config)

chr <- strsplit(config["chromosomes"], " ", fixed=TRUE)[[1]]
files <- sapply(chr, function(c) insertChromString(config["corr_file"], c, "corr_file"))

corr <- do.call(rbind, lapply(unname(files), function(f) {
    c <- getobj(f)
    dat <- data.frame(t(c$snpcorr))
    n_pcs <- min(as.integer(config["n_pcs"]), ncol(dat))
    dat <- dat[,1:n_pcs]
    names(dat) <- 1:n_pcs
    dat <- cbind(dat, chr=c$chromosome, pos=c$position, stringsAsFactors=FALSE)
}))
n_pcs <- ncol(corr) - 2

corr <- gather(corr, PC, value, -chr, -pos) %>%
    filter(!is.na(value)) %>%
    mutate(chr=factor(chr, levels=c(1:22, "X")),
           value=abs(value))

chr <- levels(corr$chr)
cmap <- setNames(rep_len(brewer.pal(8, "Dark2"), length(chr)), chr)

# plot over multiple pages
n_plots <- ceiling(n_pcs/as.integer(config["n_perpage"]))
bins <- as.integer(cut(1:n_pcs, n_plots))
for (i in 1:n_plots) {
    bin <- which(bins == i)
    p <- filter(corr, PC %in% bin) %>%
        mutate(PC=factor(paste0("PC", PC), levels=paste0("PC", 1:n_pcs))) %>%
        ggplot(aes(chr, value, group=interaction(chr, pos), color=chr)) +
        geom_point(position=position_dodge(0.8)) +
        facet_wrap(~PC, scales="free", ncol=1) +
        scale_color_manual(values=cmap, breaks=names(cmap)) +
        theme_bw() +
        theme(legend.position="none") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        xlab("chromosome") + ylab("abs(correlation)")
    ggsave(paste0(config["out_prefix"], "_" , i, ".png"), plot=p, width=10, height=15)
}
