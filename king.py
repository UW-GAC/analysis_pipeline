#! /usr/bin/env python2.7

"""Identity By Descent"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
Identity by Descent with the following steps:
1) Convert GDS to BED for use with KING
2) KING --ibdseg to get initial kinship estimates
3) KING --kinship to get divergence matrix
"""

parser = ArgumentParser(description=description)
parser.add_argument("config_file", help="configuration file")
parser.add_argument("--cluster_type", default="UW_Cluster",
                    help="type of compute cluster environment [default %(default)s]")
parser.add_argument("--cluster_file", default=None,
                    help="json file containing options to pass to the cluster")
parser.add_argument("--verbose", action="store_true", default=False,
                    help="enable verbose output to help debug")
parser.add_argument("-n", "--ncores", default="8",
                    help="number of cores to use [default %(default)s]")
parser.add_argument("-e", "--email", default=None,
                    help="email address for job reporting")
parser.add_argument("--print_only", action="store_true", default=False,
                    help="print qsub commands without submitting")
parser.add_argument("--version", action="version",
                    version="TopmedPipeline "+TopmedPipeline.__version__,
                    help="show the version number and exit")
args = parser.parse_args()

configfile = args.config_file
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

# include all samples in BED
config = deepcopy(configdict)
config["sample_include_file"] = "NA"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)
jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], email=email, print_only=print_only)


## this swaps alleles to make KING --ibdseg happy about reading the file
job = "plink_make-bed"

bedprefix = configdict["bed_file"]
arglist = ["--bfile", bedprefix, "--make-bed", "--out", bedprefix]

plinkid = cluster.submitJob(binary=True, job_name=job, cmd="plink", args=arglist, holdid=[jobid], email=email, print_only=print_only)

## when input and output files have same name, plink renames input with "~"
tmpid = cluster.submitJob(binary=True, job_name="rm_tmp_files", cmd="rm", args=[bedprefix + ".*~"], holdid=[plinkid], email=email, print_only=print_only)



job = "king_ibdseg"

bedfile = configdict["bed_file"] + ".bed"
outprefix = configdict["data_prefix"] + "_king_ibdseg"
arglist = ["-b", bedfile, "--cpus", ncores, "--ibdseg", "--prefix", outprefix]

segid = cluster.submitJob(binary=True, job_name=job, cmd="king", args=arglist, holdid=[plinkid], request_cores=ncores, email=email, print_only=print_only)

kingfile = outprefix + ".seg"


# gzip output
segid = cluster.submitJob(binary=True, job_name="gzip", cmd="gzip", args=[kingfile], holdid=[segid], email=email, print_only=print_only)
kingfile = kingfile + ".gz"


job = "kinship_plots"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["kinship_file"] = kingfile
config["kinship_method"] = "king_ibdseg"
config["out_file_all"] = configdict["plots_prefix"] + "_king_ibdseg_kinship_all.pdf"
config["out_file_cross"] = configdict["plots_prefix"] + "_king_ibdseg_kinship_cross.pdf"
config["out_file_study"] = configdict["plots_prefix"] + "_king_ibdseg_kinship_study.pdf"
configfile = configdict["config_prefix"] + "_" + job + "_ibdseg.config"
TopmedPipeline.writeConfig(config, configfile)

segplotid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[segid], email=email, print_only=print_only)


job = "king_to_matrix"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["king_file"] = kingfile
config["kinship_method"] = "king_ibdseg"
config["out_prefix"] = configdict["data_prefix"] + "_king_ibdseg_Matrix"
configfile = configdict["config_prefix"] + "_" + job + "_ibdseg.config"
TopmedPipeline.writeConfig(config, configfile)

segmatid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[segid], email=email, print_only=print_only)



job = "ibd_king"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_file"] = configdict["data_prefix"] + "_king_robust.gds"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

kinid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[jobid], request_cores=ncores, email=email, print_only=print_only)


## fails for large sample sizes (~50k)
## https://github.com/zhengxwen/SNPRelate/issues/57

# job = "kinship_plots"

# rscript = os.path.join(pipeline, "R", job + ".R")

# config = deepcopy(configdict)
# config["kinship_file"] = configdict["data_prefix"] + "_king_robust.gds"
# config["kinship_method"] = "king"
# config["out_file_all"] = configdict["plots_prefix"] + "_king_robust_kinship_all.pdf"
# config["out_file_study"] = configdict["plots_prefix"] + "_king_robust_kinship_study.pdf"
# configfile = configdict["config_prefix"] + "_" + job + "_robust.config"
# TopmedPipeline.writeConfig(config, configfile)

# kinid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[kinid], email=email, print_only=print_only)



# post analysis
job = "post_analysis"
jobpy = job + ".py"
pcmd=os.path.join(pipeline, jobpy)
argList = [pcmd, "-a", cluster.getAnalysisName(), "-l", cluster.getAnalysisLog(),
           "-s", cluster.getAnalysisStartSec()]
pdriver=os.path.join(pipeline, "run_python.sh")
holdlist = [segplotid, segmatid, kinid]
cluster.submitJob(job_name=job, cmd=pdriver, args=argList, holdid=holdlist, print_only=print_only)
