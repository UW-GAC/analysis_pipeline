# script to install R packages needed for pipeline
# must be run from analysis_pipeline directory
# argument to script is path to R library for install

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) .libPaths(args[1])

source("https://bioconductor.org/biocLite.R")
biocLite(c("SeqVarTools", "SNPRelate", "GENESIS", "survey", "CompQuadForm",
           "argparser", "ggplot2", "GGally", "rmarkdown", "devtools"),
         ask=FALSE)

devtools::install_git("git://github.com/zhengxwen/gdsfmt", dependencies=FALSE)
devtools::install_git("git://github.com/zhengxwen/SeqArray", dependencies=FALSE)
devtools::install_git("git://github.com/zhengxwen/SNPRelate", dependencies=FALSE)
devtools::install_git("git://github.com/smgogarten/GWASTools", dependencies=FALSE)
devtools::install_git("git://github.com/smgogarten/SeqVarTools", dependencies=FALSE)
devtools::install_git("git://github.com/smgogarten/GENESIS", dependencies=FALSE)

devtools::install("TopmedPipeline", dependencies=FALSE)
