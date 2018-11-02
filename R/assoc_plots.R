library(argparser)
library(TopmedPipeline)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
sessionInfo()

argp <- arg_parser("Association plots")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("assoc_file",
              "assoc_type")
optional <- c("chromosomes"="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X",
              "known_hits_file"=NA,
              "out_file_manh"="manhattan.png",
              "out_file_qq"="qq.png",
              "thin"=TRUE)
config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(basename(argv$config), ".assoc_plots.params"))

chr <- strsplit(config["chromosomes"], " ", fixed=TRUE)[[1]]
files <- sapply(chr, function(c) insertChromString(config["assoc_file"], c, "assoc_file"))

assoc <- getAssoc(files, config["assoc_type"])

## omit known hits?
if (!is.na(config["known_hits_file"]) & config["assoc_type"] == "single") {
    hits <- getobj(config["known_hits_file"])
    assoc <- omitKnownHits(assoc, hits, flank=500)
}

if ("stat" %in% names(assoc)) {
    ## burden or single
    lambda <- calculateLambda((assoc$stat)^2, df=1)
} else {
    ## SKAT
    lambda <- calculateLambda(qchisq(assoc$pval, df=1, lower=FALSE), df=1)
}

## qq plot
n <- nrow(assoc)
x <- 1:n
dat <- data.frame(obs=sort(assoc$pval),
                  exp=x/n,
                  upper=qbeta(0.025, x, rev(x)),
                  lower=qbeta(0.975, x, rev(x)))

if (as.logical(config["thin"])) {
    dat <- mutate(dat, logp=-log10(obs)) %>%
        thinPoints("logp", n=10000, nbins=10)
}

p <- ggplot(dat, aes(-log10(exp), -log10(obs))) +
    geom_line(aes(-log10(exp), -log10(upper)), color="gray") +
    geom_line(aes(-log10(exp), -log10(lower)), color="gray") +
    geom_point() +
    geom_abline(intercept=0, slope=1, color="red") +
    xlab(expression(paste(-log[10], "(expected P)"))) +
    ylab(expression(paste(-log[10], "(observed P)"))) +
    ggtitle(paste("lambda =", format(lambda, digits=4, nsmall=3))) +
    theme_bw() +
    theme(plot.title = element_text(size = 22))
ggsave(config["out_file_qq"], plot=p, width=6, height=6)

rm(dat)


## manhattan plot
chr <- levels(assoc$chr)
cmap <- setNames(rep_len(brewer.pal(8, "Dark2"), length(chr)), chr)

# significance level
if (config["assoc_type"] == "single") {
    ## genome-wide significance
    signif <- c(5e-8, 5e-9, 1e-9)
} else {
    ## bonferroni
    signif <- 0.05/nrow(assoc)
}

if (as.logical(config["thin"])) {
    assoc <- mutate(assoc, logp=-log10(pval)) %>%
        thinPoints("logp", n=10000, nbins=10, groupBy="chr")
}

p <- ggplot(assoc, aes(chr, -log10(pval), group=interaction(chr, pos), color=chr)) +
    geom_point(position=position_dodge(0.8)) +
    scale_color_manual(values=cmap, breaks=names(cmap)) +
    geom_hline(yintercept=-log10(signif), linetype='dashed') +
    theme_bw() +
    theme(legend.position="none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("chromosome") + ylab(expression(-log[10](p)))
ggsave(config["out_file_manh"], plot=p, width=10, height=5)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
