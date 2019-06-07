library(argparser)
library(TopmedPipeline)
library(rtracklayer)
sessionInfo()

argp <- arg_parser("liftover bim file")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("bim_file",
              "chain_file")
optional <- c("out_bim_file"="liftover.bim")
config <- setConfigDefaults(config, required, optional)
print(config)

bim <- read.table(config["bim_file"], as.is=TRUE)
gr <- GRanges(seqnames=paste0("chr",bim[,1]), ranges=IRanges(start=bim[,4], end=bim[,4]))

chain <- import.chain(config["chain_file"])
lift <- liftOver(gr, chain)

# report out number of successful liftovers
one.rslt <- sum(elementNROWS(lift) == 1)
no.rslt <- sum(elementNROWS(lift) == 0)
message(one.rslt, " input positions successfully converted; ", no.rslt, " positions failed conversion (will be set to NA)")
pct <- round(one.rslt/nrow(bim), 6) * 100
message(pct, "% successful conversions")

# extract converted chrom
lift.chr <- gsub("chr", "", as.character(seqnames(lift)))
lift.chr[is.na(lift.chr)] <- "0"
lift.pos <- as.integer(start((lift)))
lift.pos[is.na(lift.pos)] <- 0

bim[,1] <- lift.chr
bim[,4] <- lift.pos

write.table(bim, file=config["out_bim_file"], quote=FALSE, col.names=FALSE, row.names=FALSE)
