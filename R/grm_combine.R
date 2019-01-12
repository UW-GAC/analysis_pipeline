library(argparser)
library(TopmedPipeline)
library(SNPRelate)
sessionInfo()

argp <- arg_parser("Combine variants in files split by chromosome")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("in_file")
optional <- c("chromosomes"="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22",
              "out_file"="grm.gds")
config <- setConfigDefaults(config, required, optional)
print(config)

chr <- strsplit(config["chromosomes"], " ", fixed=TRUE)[[1]]
files <- sapply(chr, function(c) insertChromString(config["in_file"], c, "in_file"))

## write to the scratch disk of each node
outfile.tmp <- tempfile()
message("gds temporarily located at ", outfile.tmp)

snpgdsMergeGRM(files, out.fn=outfile.tmp)

## copy it
file.copy(outfile.tmp, config["out_file"])
## remove the tmp file
file.remove(outfile.tmp)

# delete by-chrom files
unlink(files)

ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
