# script to install R packages needed for pipeline
# must be run from analysis_pipeline directory
# argument to script is path to R library for install

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) .libPaths(args[1])

install.packages(c("BiocManager", "remotes"), repos="https://cloud.r-project.org")

BiocManager::install(c("SeqVarTools", "SNPRelate", "GENESIS", "survey", "CompQuadForm",
                       "argparser", "ggplot2", "GGally", "rmarkdown"),
                     ask=FALSE)
                     
remotes::install_git("git://github.com/zhengxwen/SeqArray", dependencies=FALSE)
remotes::install_git("git://github.com/smgogarten/GWASTools", dependencies=FALSE)
remotes::install_git("git://github.com/smgogarten/SeqVarTools", dependencies=FALSE)
remotes::install_git("git://github.com/UW-GAC/GENESIS", dependencies=FALSE)

remotes::install_local("TopmedPipeline", dependencies=FALSE)
