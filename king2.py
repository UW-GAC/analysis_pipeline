#! /usr/bin/env python2.7

"""Identity By Descent"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
Identity by Descent with the following steps:
1) Select SNPs with LD pruning
2) IBD calculations with KING-robust
"""

parser = ArgumentParser(description=description)
parser.add_argument("config_file", help="configuration file")
parser.add_argument("-c", "--chromosomes", default="1-22",
                    help="range of chromosomes [default %(default)s]")
parser.add_argument("--cluster_type", default="UW_Cluster",
                    help="type of compute cluster environment [default %(default)s]")
parser.add_argument("--cluster_file", default=None,
                    help="json file containing options to pass to the cluster")
parser.add_argument("--verbose", action="store_true", default=False,
                    help="enable verbose output to help debug")
parser.add_argument("-n", "--ncores", default="1-8",
                    help="number of cores to use; either a number (e.g, 1) or a range of numbers (e.g., 1-4) [default %(default)s]")
parser.add_argument("-e", "--email", default=None,
                    help="email address for job reporting")
parser.add_argument("--print_only", action="store_true", default=False,
                    help="print qsub commands without submitting")
parser.add_argument("--version", action="version",
                    version="TopmedPipeline "+TopmedPipeline.__version__,
                    help="show the version number and exit")
args = parser.parse_args()

configfile = args.config_file
chromosomes = args.chromosomes
cluster_file = args.cluster_file
cluster_type = args.cluster_type
ncores = args.ncores
email = args.email
print_only = args.print_only
verbose = args.verbose

version = "--version " + TopmedPipeline.__version__

cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type, cluster_file, verbose)

pipeline = cluster.getPipelinePath()
driver = os.path.join(pipeline, "runRscript.sh")

configdict = TopmedPipeline.readConfig(configfile)
configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "data", "log", "plots"])

# analysis init
cluster.analysisInit(print_only=print_only)


job = "gds2bed"

rscript = os.path.join(pipeline, "R", job + ".R")

jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], email=email, print_only=print_only)


## this swaps alleles to make KING --ibdseg happy about reading the file
job = "plink_make-bed"

bedprefix = configdict["bed_file"]
#arglist = ["--bfile", bedprefix, "--make-bed", "--out", "data/tmp"]
arglist = ["--bfile", bedprefix, "--make-bed", "--out", bedprefix]

job_cmd = cluster.clusterCfg["submit_cmd"]
subOpts = deepcopy(cluster.clusterCfg["submit_opts"])
subOpts["-b"] = "y"
jobid = cluster.executeJobCmd(subOpts, job_cmd=job_cmd, job_name=job, cmd="plink", args=arglist, holdid=[jobid], email=email, print_only=print_only)
tmpid = cluster.executeJobCmd(subOpts, job_cmd=job_cmd, job_name="rm_tmp_files", cmd="rm", args=[bedprefix + ".*~"], holdid=[jobid], email=email, print_only=print_only)

## overwrite original BED
# job = "rename_files"
# outF = open("log/tmp.sh", "w")
# outF.write("#! /bin/bash\n")
# for ext in [".bed", ".bim", ".fam"]:
#     outF.write("mv tmp" + ext + " " + bedprefix + ext + "\n")
# outF.close()
# jobid = cluster.submitJob(job_name=job, cmd="log/tmp.sh", email=email, print_only=print_only)


job = "king_ibdseg"

bedfile = configdict["bed_file"] + ".bed"
outprefix = configdict["data_prefix"] + "_king"
arglist = ["-b", bedfile, "--lessmem", "--ibdseg", "--prefix", outprefix]

# job_cmd = cluster.clusterCfg["submit_cmd"]
# subOpts = deepcopy(cluster.clusterCfg["submit_opts"])
# subOpts["-b"] = "y"
jobid = cluster.executeJobCmd(subOpts, job_cmd=job_cmd, job_name=job, cmd="king", args=arglist, holdid=[jobid], email=email, print_only=print_only)


## format ibdseg results with kingToMatrix
## set sample.include to all samples so every sample is represented in output sparse matrix
## set threshold to 2^(-13/2) (5th deg)


## run ibd_king R script on GDS file

## modify kinship_plots to use hexbin instead of geom_point
## option to read text file with king ibdseg output



# job = "kinship_plots"

# rscript = os.path.join(pipeline, "R", job + ".R")

# config = deepcopy(configdict)
# config["kinship_file"] = configdict["data_prefix"] + "_ibd_king.gds"
# config["kinship_method"] = "king"
# config["out_file_all"] = configdict["plots_prefix"] + "_kinship_all.pdf"
# config["out_file_cross"] = configdict["plots_prefix"] + "_kinship_cross.pdf"
# config["out_file_study"] = configdict["plots_prefix"] + "_kinship_study.pdf"
# configfile = configdict["config_prefix"] + "_" + job + ".config"
# TopmedPipeline.writeConfig(config, configfile)

# jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[jobid], email=email, print_only=print_only)


# post analysis
job = "post_analysis"
jobpy = job + ".py"
pcmd=os.path.join(pipeline, jobpy)
argList = [pcmd, "-a", cluster.getAnalysisName(), "-l", cluster.getAnalysisLog(),
           "-s", cluster.getAnalysisStartSec()]
pdriver=os.path.join(pipeline, "run_python.sh")
cluster.submitJob(job_name=job, cmd=pdriver, args=argList,
                  holdid=[jobid], print_only=print_only)
