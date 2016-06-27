library(argparser)
library(TopmedPipeline)
library(Biobase)
library(GENESIS)
library(gdsfmt)
sessionInfo()

argp <- arg_parser("Null model for association tests")
argp <- add_argument(argp, "config", help="path to config file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("outcome",
              "pca_file",
              "pcrelate_file",
              "phenotype_file")
optional <- c("binary"=FALSE,
              "covars"=NA,
              "inverse_normal"=FALSE,
              "n_pcs"=3,
              "out_file"="null_model.RData",
              "sample_include_file"=NA)
config <- setConfigDefaults(config, required, optional)
print(config)

# get PCs
pca <- getobj(config["pca_file"])
n_pcs <- as.integer(config["n_pcs"])
pcs <- pca$vectors[,1:n_pcs]
colnames(pcs) <- paste0("PC", 1:n_pcs)

# select samples
if (!is.na(config["sample_include_file"])) {
    sample.id <- getobj(config["sample_include_file"])
} else {
    sample.id <- rownames(pcs)
}

# get phenotypes
annot <- getobj(config["phenotype_file"])
annot <- annot[annot$sample.id %in% sample.id,]

# use ordering of phenotype file (should match GDS)
sample.id <- annot$sample.id
dat <- cbind(pData(annot), pcs[as.character(sample.id),])
pData(annot) <- dat

# load GRM for selected samples only
pcr <- openfn.gds(config["pcrelate_file"])
grm <- pcrelateMakeGRM(pcr, scan.include=sample.id)
closefn.gds(pcr)

# outcome and covariates
outcome <- config["outcome"]
if (!is.na(config["covars"])) {
    covars <- strsplit(config["covars"], " ", fixed=TRUE)[[1]]
} else {
    covars <- NULL
}
covars <- c(covars, colnames(pcs))

if (as.logical(config["binary"])) {
    stopifnot(all(dat[[outcome]] %in% c(0,1)))
    family <- binomial
} else {
    family <- gaussian
}

message("Model: ", outcome, " ~ ", paste(c(covars, "(1|kinship)"), collapse=" + "))
nullmod <- fitNullMM(annot, outcome=outcome, covars=covars,
                     covMatList=grm, scan.include=sample.id,
                     family=family)

# if we need an inverse normal transform, take residuals and refit null model
# for second model fit, use PCs and kinship but not other covariates
if (as.logical(config["inverse_normal"])) {
    resid.norm <- rankNorm(nullmod$resid.marginal)
    annot$resid.norm <- resid.norm[match(annot$sample.id, nullmod$scanID)]
    message(paste0("resid.norm = rankNorm(resid.marginal(", outcome, " ~ ", paste(c(covars, "(1|kinship)"), collapse=" + "), "))"))
    message("Model: resid.norm ~ ", paste(c(colnames(pcs), "(1|kinship)"), collapse=" + "))
    nullmod <- fitNullMM(annot, outcome="resid.norm", covars=colnames(pcs),
                         covMatList=grm, scan.include=sample.id,
                         family=family)
}

save(nullmod, file=config["out_file"])
