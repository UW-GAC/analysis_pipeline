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

devtools::install_git("git://github.com/r-lib/devtools.git")
library(git2r)
library(devtools)
install_git("git://github.com/smgogarten/GWASTools.git", branch="v1.27.1", dependencies=FALSE)
install_git("git://github.com/smgogarten/SeqVarTools.git", branch="v1.19.2", dependencies=FALSE)
install_git("git://github.com/smgogarten/GENESIS.git", branch="v2.11.5", dependencies=FALSE)

install("TopmedPipeline", dependencies=FALSE)
