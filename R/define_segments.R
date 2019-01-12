library(argparser)
library(TopmedPipeline)
sessionInfo()

argp <- arg_parser("Define genome segments for association test")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--segment_length", help="segment length in kb", type="integer", default=10000)
argp <- add_argument(argp, "--n_segments", help="number of segments (overrides segment length)", type="integer", default=NA)
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
seg.length <- argv$segment_length * 1000
n <- argv$n_segments

required <- c()
optional <- c("genome_build"="hg38",
              "out_file"="segments.txt")
config <- setConfigDefaults(config, required, optional)
print(config)

if (!is.na(n)) {
    segments <- defineSegments(n=n, build=config["genome_build"])
} else {
    segments <- defineSegments(seg.length=seg.length, build=config["genome_build"])
}

writeSegmentFile(segments, config["out_file"])
