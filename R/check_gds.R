library(argparser)
library(TopmedPipeline)
library(SeqArray)
library(tools)
sessionInfo()

argp <- arg_parser("Check GDS against original VCF")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)

required <- c("vcf_file", "gds_file")
optional <- c("sample_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

## vcf file can have two parts split by chromosome identifier
vcffile <- config["vcf_file"]
gdsfile <- config["gds_file"]
if (!is.na(chr)) {
    vcffile <- insertChromString(vcffile, chr, "vcf_file")
    gdsfile <- insertChromString(gdsfile, chr, "gds_file")
}

# compare genotypes via pipes

# no INFO and FORMAT, except genotypes
cmd_gds2vcf <- "
suppressPackageStartupMessages(library(SeqArray))
f <- seqOpen('%s')
on.exit(seqClose(f))
%s
seqGDS2VCF(f, stdout(), info.var=character(), fmt.var=character(), verbose=FALSE)
"
if (!is.na(config["sample_file"])) {
    keep_str <- sprintf("
sample.id <- readLines('%s')
seqSetFilter(f, sample.id=sample.id)
", config["sample_file"])
} else {
    keep_str <- ""
}
tmpRscript <- tempfile()
writeLines(sprintf(cmd_gds2vcf, gdsfile, keep_str), tmpRscript)
cmd_gds <- sprintf("Rscript --vanilla %s", tmpRscript)
f1 <- pipe(cmd_gds, "rt")
message("    open: ", basename(gdsfile), "\n")

cmd_bcftools <- "bcftools annotate -x INFO,^FORMAT/GT %s"
cmd_vcf <- sprintf(cmd_bcftools, vcffile)
f2 <- pipe(cmd_vcf, "rt")
message("    open: ", basename(vcffile), "\n")


# skip header
l1 <- 1L
while (length(s <- readLines(f1, n=1L))) {
    l1 <- l1 + 1L
    if (substr(s, 1L, 6L) == "#CHROM") break
}
message(sprintf("GDS to VCF text, start from line %d\n", l1))

l2 <- 1L
while (length(s <- readLines(f2, n=1L))) {
    l2 <- l2 + 1L
    if (substr(s, 1L, 6L) == "#CHROM") break
}
message(sprintf("BCF to VCF text, start from line %d\n", l2))

l1 <- l2 <- 0L
iline <- 0L

nmax <- 10000L

# for-loop
repeat {
    s1 <- readLines(f1, n=nmax)
    s2 <- readLines(f2, n=nmax)
    if (length(s2) > 0L & substr(s2[1], 1, 3) == "chr") {
        s2 <- sub("^chr", "", s2) # chr2 in BCF, but 2 in GDS, strip off prefix "chr"
    }
    if (!identical(s1, s2)) {
        if (length(s1) < length(s2)) {
            i <- which(s1 != s2[seq_len(length(s1))])[1]
            if (!is.na(i)) {
                message("GDS line:\n", s1[i], "\n")
                message("BCF line:\n", s2[i], "\n")
                stop("GDS line ", l1+i, " is different from BCF line ", l2+i, ".")
            } else {
                message("BCF line:\n", s2[length(s1)+1], "\n")
                stop("GDS has ", l1+length(s1), " lines, but BCF has more.")
            }
        } else if (length(s1) > length(s2)) {
            i <- which(s1[seq_len(length(s2))] != s2)[1]
            if (!is.na(i)) {
                message("GDS line:\n", s1[i], "\n")
                message("BCF line:\n", s2[i], "\n")
                stop("GDS line ", l1+i, " is different from BCF line ", l2+i, ".")
            } else {
                message("GDS line:\n", s1[length(s2)+1], "\n")
                stop("BCF has ", l2+length(s2), " lines, but GDS has more.")
            }
        }
        i <- which(s1 != s2)[1]
        message("GDS line:\n", s1[i], "\n")
        message("BCF line:\n", s2[i], "\n")
        stop("GDS line ", l1+i, " is different from BCF line ", l2+i, ".")
    }
    if (length(s1)==0L | length(s2)==0L)
        break
    l1 <- l1 + length(s1)
    l2 <- l2 + length(s2)
    iline <- iline + 1L
    if ((iline >= 100L) | (length(s1) < nmax)) {
        message("    passed: ", l1, " lines\n")
        iline <- 0L
    }
}

close(f1)
close(f2)
unlink(tmpRscript)
