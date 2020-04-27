#! /usr/bin/env python3

"""Identity By Descent"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
Identity by Descent with the following steps:
1) KING robust using SNPRelate algorithm
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
submitPath = cluster.getSubmitPath()
driver = os.path.join(submitPath, "runRscript.sh")

configdict = TopmedPipeline.readConfig(configfile)
configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "data", "log", "plots"])

# analysis init
cluster.analysisInit(print_only=print_only)


job = "ibd_king"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_file"] = configdict["data_prefix"] + "_king_robust.gds"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

kinid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], request_cores=ncores, email=email, print_only=print_only)


job = "kinship_plots"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["kinship_file"] = configdict["data_prefix"] + "_king_robust.gds"
config["kinship_method"] = "king"
config["out_file_all"] = configdict["plots_prefix"] + "_king_robust_kinship_all.pdf"
config["out_file_study"] = configdict["plots_prefix"] + "_king_robust_kinship_study.pdf"
configfile = configdict["config_prefix"] + "_" + job + "_robust.config"
TopmedPipeline.writeConfig(config, configfile)

kinid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[kinid], email=email, print_only=print_only)


# post analysis
bname = "post_analysis"
job = "king" + "_" + bname
jobpy = bname + ".py"
pcmd=os.path.join(submitPath, jobpy)
holdlist = [kinid]
argList = ["-a", cluster.getAnalysisName(), "-l", cluster.getAnalysisLog(),
           "-s", cluster.getAnalysisStartSec()]
cluster.submitJob(binary=True, job_name=job, cmd=pcmd, args=argList,
                  holdid=holdlist, print_only=print_only)
