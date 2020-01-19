#! /usr/bin/env python3

"""PC-Relate"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
PC-Relate with the following steps:
1) Calculate ISAF betas
2) Calculate kinship in all sample blocks in parallel
3) Correct kinship estimates
"""

parser = ArgumentParser(description=description)
parser.add_argument("config_file", help="configuration file")
parser.add_argument("--cluster_type", default="UW_Cluster",
                    help="type of compute cluster environment [default %(default)s]")
parser.add_argument("--cluster_file", default=None,
                    help="json file containing options to pass to the cluster")
parser.add_argument("--verbose", action="store_true", default=False,
                    help="enable verbose output to help debug")
parser.add_argument("-e", "--email", default=None,
                    help="email address for job reporting")
parser.add_argument("--print_only", action="store_true", default=False,
                    help="print cluster commands without submitting")
parser.add_argument("--version", action="version",
                    version="TopmedPipeline "+TopmedPipeline.__version__,
                    help="show the version number and exit")
args = parser.parse_args()

configfile = args.config_file
cluster_file = args.cluster_file
cluster_type = args.cluster_type
email = args.email
print_only = args.print_only
verbose = args.verbose

version = "--version " + TopmedPipeline.__version__

cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type, cluster_file, verbose)

pipeline = cluster.getPipelinePath()
submitPath = cluster.getSubmitPath()
driver = os.path.join(submitPath, "runRscript.sh")

configdict = TopmedPipeline.readConfig(configfile)
configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "data", "log", "plots"])

# analysis init
cluster.analysisInit(print_only=print_only)

job = "pcrelate_beta"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_prefix"] = configdict["data_prefix"] + "_isaf_beta"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], email=email, print_only=print_only)


job = "pcrelate"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["beta_file"] = configdict["data_prefix"] + "_isaf_beta.RData"
config["out_prefix"] = configdict["data_prefix"]
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

# define sample blocks
n_blocks = int(configdict.setdefault("n_sample_blocks", "1"))
if n_blocks > 1:
    blocks = [(i,j) for i in range(1, n_blocks+1) for j in range(i, n_blocks+1)]
    block_range = "1-" + str(len(blocks))
else:
    block_range = "1"

jobid = cluster.submitJob(job_name=job, cmd=driver, args=["-s", rscript, configfile, version], holdid=[jobid], array_range=block_range, email=email, print_only=print_only)


job = "pcrelate_correct"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["pcrelate_prefix"] = configdict["data_prefix"]
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[jobid], email=email, print_only=print_only)


job = "kinship_plots"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["kinship_file"] = configdict["data_prefix"] + "_pcrelate.RData"
config["kinship_method"] = "pcrelate"
config["out_file_all"] = configdict["plots_prefix"] + "_kinship_all.pdf"
config["out_file_cross"] = configdict["plots_prefix"] + "_kinship_cross.pdf"
config["out_file_study"] = configdict["plots_prefix"] + "_kinship_study.pdf"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[jobid], email=email, print_only=print_only)

# post analysis
bname = "post_analysis"
job = "pcrelate" + "_" + bname
jobpy = bname + ".py"
pcmd=os.path.join(submitPath, jobpy)
argList = ["-a", cluster.getAnalysisName(), "-l", cluster.getAnalysisLog(),
           "-s", cluster.getAnalysisStartSec()]
cluster.submitJob(binary=True, job_name=job, cmd=pcmd, args=argList,
                  holdid=[jobid], print_only=print_only)
