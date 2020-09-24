# script to install R packages needed for pipeline
# must be run from analysis_pipeline directory
# argument to script is path to R library for install

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) .libPaths(args[1])

install.packages(c("BiocManager", "remotes"), repos="https://cloud.r-project.org")

BiocManager::install(c("SeqVarTools", "SNPRelate", "GENESIS", "survey", "CompQuadForm",
                       "argparser", "data.table", "ggplot2", "GGally", "hexbin",
                       "lazyeval", "logistf", "R.utils", "rmarkdown", "SPAtest"),
                     update=FALSE, ask=FALSE)

# if R version is current, BiocManager will automatically install latest release
if (getRversion() < "3.6.0") {
    remotes::install_git("git://github.com/zhengxwen/gdsfmt.git", ref="v1.20.0", dependencies=FALSE)
    remotes::install_git("git://github.com/zhengxwen/SeqArray.git", ref="v1.24.0", dependencies=FALSE)
    remotes::install_git("git://github.com/smgogarten/GWASTools.git", ref="v1.32.0", dependencies=FALSE)
}

remotes::install_git("git://github.com/smgogarten/SeqVarTools.git", ref="v1.27.1", dependencies=FALSE)
remotes::install_git("git://github.com/UW-GAC/GENESIS.git", ref="v2.19.5", dependencies=FALSE)

remotes::install_local("TopmedPipeline", dependencies=FALSE)

