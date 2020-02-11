library(argparser)
library(TopmedPipeline)
library(Biobase)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Partition samples into related and unrelated sets")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("kinship_file")
optional <- c("divergence_file"=NA,
              "kinship_threshold"=0.04419417, # 2^(-9/2), 3rd degree
              "divergence_threshold"=-0.04419417, # 2^(-9/2), 3rd degree
              "out_related_file"="related.RData",
              "out_unrelated_file"="unrelated.RData",
              "phenotype_file"=NA,
              "sample_include_file"=NA,
              "study"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
    message("Using ", length(sample.id), " samples")
} else {
    sample.id <- NULL
    message("Using all samples")
}

kinMat <- kinobj(config["kinship_file"])

if (!is.na(config["divergence_file"])) {
    divMat <- kinobj(config["divergence_file"])
    message("Using divergence matrix to find unrelated set")
} else {
    divMat <- NULL
    message("No divergence matrix specified")
}

kin_thresh <- as.numeric(config["kinship_threshold"])
div_thresh <- as.numeric(config["divergence_threshold"])

## for each study, find median kinship coefficient
## if median KC > 0, run find_unrelated on that study only, with threshold of 2^(-9/2) + median KC
if (!is.na(config["study"]) & !is.na(config["phenotype_file"])) {
    # could replace this with a read.gdsn command, but unlikely we will need it
    if (is(kinMat, "gds.class")) stop("can't compute median kinship by study on GDS object")
    message("Computing median kinship by study")
    
    study <- config["study"]
    annot <- getobj(config["phenotype_file"])
    stopifnot(study %in% varLabels(annot))
    annot <- pData(annot)[,c("sample.id", study)]
    names(annot)[2] <- "study"
    if (!is.null(sample.id)) {
        annot <- annot[annot$sample.id %in% sample.id,]
    }
    studies <- unique(annot$study)
    study.partition <- list()
    for (s in studies) {
        ids <- annot$sample.id[annot$study %in% s]
        ind <- rownames(kinMat) %in% ids
        medKC <- medianKinship(kinMat[ind,ind])
        if (medKC > 0) {
            message("Median kinship for ", s, " is ", medKC, ".\n",
                    "Finding unrelated set separately.")
            study.partition[[s]] <- pcairPartition(kinobj=kinMat, kin.thresh=(kin_thresh + medKC),
                                                   divobj=divMat, div.thresh=div_thresh,
                                                   sample.include=ids)
            message("Found ", length(study.partition[[s]]$unrels), " unrelated and ", length(study.partition[[s]]$rels), " related samples")
        }
    }
    
    ## combine unrelated samples from individual studies
    # will be NULL if list is empty
    study.unrel <- unlist(lapply(study.partition, function(x) x$unrels), use.names=FALSE)
} else {
    study.unrel <- NULL
}


## run pcairPartition on everyone, forcing list of per-study unrel into unrelated set
part <- pcairPartition(kinobj=kinMat, kin.thresh=kin_thresh,
                       divobj=divMat, div.thresh=div_thresh,
                       sample.include=sample.id,
                       unrel.set=study.unrel)

rels <- part$rels
unrels <- part$unrels
save(rels, file=config["out_related_file"])
save(unrels, file=config["out_unrelated_file"])
message("Found ", length(unrels), " unrelated and ", length(rels), " related samples")

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
