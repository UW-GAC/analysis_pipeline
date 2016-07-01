library(argparser)
library(TopmedPipeline)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
sessionInfo()

argp <- arg_parser("Association plots")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("assoc_file",
              "assoc_type")
optional <- c("chromosomes"="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X",
              "out_file_manh"="manhattan.png",
              "out_file_qq"="qq.png")
config <- setConfigDefaults(config, required, optional)
print(config)

chr <- strsplit(config["chromosomes"], " ", fixed=TRUE)[[1]]
files <- sapply(chr, function(c) insertChromString(config["assoc_file"], c, "assoc_file"))

assoc <- lapply(files, getobj)
names(assoc) <- NULL

if (config["assoc_type"] == "single") {
    assoc <- do.call(rbind, assoc)
} else if (config["assoc_type"] == "aggregate") {
    assoc <- do.call(rbind, lapply(assoc, function(x) {
        tmp <- x$results %>%
            add_rownames("group_id") %>%
            filter(n.site > 0)
        group.info <- do.call(rbind, lapply(tmp$group_id, function(g) {
            grp <- x$variantInfo[[g]][1, c("chr", "pos")]
            grp$group_id <- g
            grp
        }))
        left_join(tmp, group.info, by="group_id")
    }))
} else if (config["assoc_type"] == "window") {
    assoc <- do.call(rbind, lapply(assoc, function(x) {
        filter(x$results, n.site > 0, dup == 0) %>%
            mutate(pos=floor((window.start + window.stop)/2))
    }))
} else {
    stop("assoc_type should be 'single', 'aggregate' or 'window'")
}

if ("pval_0" %in% names(assoc)) {
    ## SKAT
    pval.col <- if ("pval_SKATO" %in% names(assoc)) "pval_SKATO" else "pval_0"
    assoc <- select_(assoc, "chr", "pos", pval.col) %>%
        rename_(pval=pval.col)
    lambda <- calculateLambda(qchisq(assoc$pval, df=1, lower=FALSE), df=1)
} else {
    ## burden or single
    assoc <- select(assoc, chr, pos, ends_with("stat"), ends_with("pval"))
    names(assoc)[3:4] <- c("stat", "pval")
    lambda <- calculateLambda(assoc$stat, df=1)
}
assoc <- filter(assoc, !is.na(pval)) %>%
    mutate(chr=factor(chr, levels=c(1:22, "X")))


## manhattan plot
chr <- levels(assoc$chr)
cmap <- setNames(rep_len(brewer.pal(8, "Dark2"), length(chr)), chr)

p <- ggplot(assoc, aes(chr, -log10(pval), group=interaction(chr, pos), color=chr)) +
    geom_point(position=position_dodge(0.8)) +
    scale_color_manual(values=cmap, breaks=names(cmap)) +
    geom_hline(yintercept=-log10(5e-8), linetype='dashed') +
    theme_bw() +
    theme(legend.position="none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("chromosome") + ylab(expression(-log[10](p)))
ggsave(config["out_file_manh"], plot=p, width=10, height=5)


## qq plot
n <- nrow(assoc)
x <- 1:n
dat <- data.frame(obs=sort(assoc$pval),
                  exp=x/n,
                  upper=qbeta(0.025, x, rev(x)),
                  lower=qbeta(0.975, x, rev(x)))
    
p <- ggplot(dat, aes(-log10(exp), -log10(obs))) +
    geom_line(aes(-log10(exp), -log10(upper)), color="gray") +
    geom_line(aes(-log10(exp), -log10(lower)), color="gray") +
    geom_point() +
    geom_abline(intercept=0, slope=1, color="red") +
    xlab(expression(paste(-log[10], "(expected P)"))) +
    ylab(expression(paste(-log[10], "(observed P)"))) +
    ggtitle(paste("lambda =", format(lambda, digits=4, nsmall=3))) +
    theme_bw()
ggsave(config["out_file_qq"], plot=p, width=6, height=6)
