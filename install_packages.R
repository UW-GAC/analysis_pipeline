# script to install R packages needed for pipeline
# must be run from analysis_pipeline directory
# argument to script is path to R library for install

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) .libPaths(args[1])

install.packages(c("BiocManager", "remotes"), repos="https://cloud.r-project.org")

BiocManager::install(c("SeqVarTools", "SNPRelate", "GENESIS", "survey", "CompQuadForm",
                       "argparser", "data.table", "ggplot2", "GGally", "hexbin", "R.utils",
                       "rmarkdown"),
                     ask=FALSE)

remotes::install_git("git://github.com/zhengxwen/gdsfmt.git", ref="v1.20.0", dependencies=FALSE)
remotes::install_git("git://github.com/zhengxwen/SeqArray.git", ref="v1.24.0", dependencies=FALSE)
remotes::install_git("git://github.com/smgogarten/GWASTools.git", dependencies=FALSE)
remotes::install_git("git://github.com/smgogarten/SeqVarTools.git", ref="v1.23.1", dependencies=FALSE)
remotes::install_git("git://github.com/UW-GAC/GENESIS.git", ref="v2.15.4", dependencies=FALSE)

remotes::install_local("TopmedPipeline", dependencies=FALSE)
