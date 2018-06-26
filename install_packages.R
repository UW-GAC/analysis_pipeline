# script to install R packages needed for pipeline
# must be run from analysis_pipeline directory
# argument to script is path to R library for install

args <- commandArgs(trailingOnly=TRUE)
.libPaths(args[1])

source("https://bioconductor.org/biocLite.R")
biocLite(c("SeqVarTools", "SNPRelate", "GENESIS", "argparser", "dplyr", "tidyr", "ggplot2", "GGally", "rmarkdown", "Matrix", "devtools"))

devtools:install_github("smgogarten/GWASTools", ref="v1.27.1", dependencies=FALSE)
devtools:install_github("smgogarten/SeqVarTools", ref="v1.19.2", dependencies=FALSE)
devtools:install_github("smgogarten/GENESIS", ref="v2.11.5", dependencies=FALSE)

devtools::install("TopmedPipeline", dependencies=FALSE)
