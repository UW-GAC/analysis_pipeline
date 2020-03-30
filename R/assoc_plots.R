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
              "plot_mac_threshold"=NA,
              "thin"=TRUE,
              "thin_npoints"=10000,
              "thin_nbins"=10,
              "truncate_pval_threshold" = 1e-12)

config <- setConfigDefaults(config, required, optional)
print(config)
writeConfig(config, paste0(basename(argv$config), ".assoc_plots.params"))

chr <- strsplit(config["chromosomes"], " ", fixed=TRUE)[[1]]
files <- sapply(chr, function(c) insertChromString(config["assoc_file"], c, "assoc_file"))

truncate_pval_threshold = as.numeric(config["truncate_pval_threshold"])

assoc <- getAssoc(files, config["assoc_type"])

## change p-value from 0 to very small number for extreme test statistics
if ("stat" %in% names(assoc)) {
    extreme <- abs(assoc$stat) > 38.5
    if (any(extreme)) assoc$pval[extreme] <- 5e-324
}

## filter plot by MAC?
if (!is.na(config["plot_mac_threshold"])) {
    mac.thresh <- as.integer(config["plot_mac_threshold"])
    assoc <- filter(assoc, MAC >= mac.thresh)
}

## omit known hits?
if (!is.na(config["known_hits_file"]) & config["assoc_type"] == "single") {
    hits <- getobj(config["known_hits_file"])
    assoc <- omitKnownHits(assoc, hits, flank=500)
}

if ("stat" %in% names(assoc)) {
    ## burden or single
    lambda <- calculateLambda((assoc$stat)^2, df=1)
    lambda_by_chr <- assoc %>%
      group_by(chr) %>%
      summarise(lambda = calculateLambda(stat^2, df=1))
} else {
    ## SKAT
    lambda <- calculateLambda(qchisq(assoc$pval, df=1, lower=FALSE), df=1)
    lambda_by_chr <- assoc %>%
      group_by(chr) %>%
      summarise(lambda = calculateLambda(qchisq(pval, df = 1, lower=FALSE), df = 1))
}

# Check if we should also generate truncated plots.
truncate = any(assoc$pval < truncate_pval_threshold)

## qq plot
dat <- assoc %>%
  select(
    obs = pval
  ) %>%
  arrange(obs) %>%
  mutate(
    x = row_number(),
    exp = x / n(),
    upper = qbeta(0.025, x, rev(x)),
    lower = qbeta(0.975, x, rev(x))
  ) %>%
  select(-x)

## Calculate QQ values separately for each chromosome.
dat_by_chr <- assoc %>%
  select(
    chr = chr,
    obs = pval
  ) %>%
  group_by(chr) %>%
  arrange(obs) %>%
  mutate(
    x = row_number(),
    exp = x / n(),
    upper = qbeta(0.025, x, rev(x)),
    lower = qbeta(0.975, x, rev(x))
  ) %>%
  select(-x) %>%
  ungroup()


thin.n <- as.integer(config["thin_npoints"])
thin.bins <- as.integer(config["thin_nbins"])
if (as.logical(config["thin"])) {
    dat <- mutate(dat, logp=-log10(obs)) %>%
        thinPoints("logp", n=thin.n, nbins=thin.bins)

    dat_by_chr <- dat_by_chr %>%
      mutate(logp = -log10(obs)) %>%
      group_by(chr) %>%
      thinPoints("logp", n=thin.n, nbins=thin.bins) %>%
      ungroup()
}

gg_qq <- list(
  geom_abline(intercept=0, slope=1, color="red"),
  geom_line(aes(-log10(exp), -log10(upper)), color="gray"),
  geom_line(aes(-log10(exp), -log10(lower)), color="gray"),
  xlab(expression(paste(-log[10], "(expected P)"))),
  ylab(expression(paste(-log[10], "(observed P)"))),
  theme_bw()
)

p <- ggplot(dat, aes(-log10(exp), -log10(obs))) +
    ggtitle(paste("lambda =", format(lambda, digits=4, nsmall=3))) +
    gg_qq +
    theme(plot.title = element_text(size = 22)) +
    geom_point()
ggsave(config["out_file_qq"], plot=p, width=6, height=6)

# QQ plots by chromosome.
p_by_chr <- ggplot(dat_by_chr, aes(-log10(exp), -log10(obs))) +
    gg_qq +
    ggtitle(paste("lambda =", format(lambda, digits=4, nsmall=3))) +
    theme(plot.title = element_text(size = 22)) +
    geom_point(size = 0.5) +
    facet_wrap(~ chr) +
    # Add lambda by chromosome.
    geom_text(data = lambda_by_chr, aes(label = sprintf("lambda == %4.3f", lambda)),
              x = -Inf, y = Inf, hjust = -0.2, vjust = 1.2, parse=T)
outfile <- gsub(".", "_bychr.", config["out_file_qq"], fixed=TRUE)
ggsave(outfile, plot = p_by_chr, width = 10, height = 9)

if (truncate) {
  p <- p +
    scale_y_continuous(
      oob = scales::squish,
      limits = c(0, -log10(truncate_pval_threshold)),
      expand = c(0,0)
    ) +
    ylab(expression(paste(-log[10], "(observed P)") ~ " (truncated)"))
  outfile <- gsub(".", "_truncated.", config["out_file_qq"], fixed=TRUE)
  ggsave(outfile, plot=p, width=6, height=6)

  # QQ plots by chromosome, truncated.
  p_by_chr <- p_by_chr +
    scale_y_continuous(
      oob = scales::squish,
      limits = c(0, -log10(truncate_pval_threshold)),
      expand = c(0,0)
    ) +
    ylab(expression(paste(-log[10], "(observed P)") ~ " (truncated)"))
  outfile <- gsub(".", "_truncated_bychr.", config["out_file_qq"], fixed=TRUE)
  ggsave(outfile, plot = p_by_chr, width = 6, height = 6)
}

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
        thinPoints("logp", n=thin.n, nbins=thin.bins, groupBy="chr")
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

if (truncate) {
  p <- p +
    scale_y_continuous(
      oob = scales::squish,
      limits = c(0, -log10(truncate_pval_threshold)),
      expand = c(0,0)
    ) +
    ylab(expression(-log[10](p) ~ " (truncated)"))
  outfile <- gsub(".", "_truncated.", config["out_file_manh"], fixed=TRUE)
  ggsave(outfile, plot=p, width=10, height=5)
}

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
