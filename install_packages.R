# script to install R packages needed for pipeline
# must be run from analysis_pipeline directory
# argument to script is path to R library for install

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) .libPaths(args[1])

install.packages(c("BiocManager", "remotes"), repos="https://cloud.r-project.org")

BiocManager::install(c("SeqVarTools", "SNPRelate", "GENESIS", "survey", "CompQuadForm",
                       "argparser", "ggplot2", "GGally", "rmarkdown"),
                     ask=FALSE)


#remotes::install_git("git://github.com/zhengxwen/SeqArray.git", dependencies=FALSE)
#remotes::install_github("zhengxwen/SeqArray",  ref="31b9556", dependencies=FALSE)
download.file("https://github.com/zhengxwen/SeqArray/archive/31b955666380fe1f92b5e2f2f56789a358c2cb77.zip",
              destfile="SeqArray.zip", method="libcurl")
system("unzip SeqArray.zip")
install.packages("SeqArray-31b955666380fe1f92b5e2f2f56789a358c2cb77", repos=NULL, type="source")
remotes::install_git("git://github.com/smgogarten/GWASTools.git", dependencies=FALSE)
remotes::install_git("git://github.com/smgogarten/SeqVarTools.git", dependencies=FALSE)
remotes::install_git("git://github.com/UW-GAC/GENESIS.git", dependencies=FALSE)

remotes::install_local("TopmedPipeline", dependencies=FALSE)
