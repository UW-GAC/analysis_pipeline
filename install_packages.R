# script to install R packages needed for pipeline
# must be run from analysis_pipeline directory
# argument to script is path to R library for install

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) .libPaths(args[1])

source("https://bioconductor.org/biocLite.R")
biocLite(c("SeqVarTools", "SNPRelate", "GENESIS",
           "survey", "CompQuadForm",
           "dplyr", "tidyr", "ggplot2", "GGally",
           "argparser", "rmarkdown", "devtools"),
         ask=FALSE)

devtools::install_url("https://github.com/smgogarten/GWASTools/archive/v1.27.1.tar.gz")
devtools::install_url("https://github.com/smgogarten/SeqVarTools/archive/v1.19.2.tar.gz")
devtools::install_url("https://github.com/smgogarten/GENESIS/archive/v2.11.5.tar.gz")

devtools::install("TopmedPipeline", dependencies=FALSE)
