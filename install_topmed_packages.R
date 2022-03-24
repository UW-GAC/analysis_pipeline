# script to install only GENESIS and TopmedPipeline
# use with caution: gets latest master branch on GitHub, not tagged release
# must be run from analysis_pipeline directory
# argument to script is path to R library for install

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) .libPaths(args[1])

remotes::install_github("UW-GAC/GENESIS", dependencies=FALSE)

remotes::install_local("TopmedPipeline", dependencies=FALSE)
