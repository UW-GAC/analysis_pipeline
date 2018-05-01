# script to install R packages needed for pipeline
# must be run from analysis_pipeline directory
# argument to script is path to R library for install

args <- commandArgs(trailingOnly=TRUE)
.libPaths(args[1])

source("https://bioconductor.org/biocLite.R")
biocLite(c("SeqVarTools", "SNPRelate", "GENESIS", "argparser", "dplyr", "tidyr", "ggplot2", "GGally", "rmarkdown", "Matrix", "devtools"))

devtools::install_github("zhengxwen/gdsfmt", dependencies=FALSE)
devtools::install_github("zhengxwen/SeqArray", dependencies=FALSE)
devtools::install_github("zhengxwen/SNPRelate", dependencies=FALSE)
devtools::install_github("smgogarten/SeqVarTools", dependencies=FALSE)
devtools::install_github("smgogarten/GENESIS", dependencies=FALSE)

devtools::install("TopmedPipeline", dependencies=FALSE)
