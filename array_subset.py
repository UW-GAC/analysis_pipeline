#! /usr/bin/env python3

"""Subset sequencing GDS files for array concordance"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
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
email = args.email
print_only = args.print_only
verbose = args.verbose

version = "--version " + TopmedPipeline.__version__

cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type, cluster_file, verbose)

pipeline = cluster.getPipelinePath()
driver = os.path.join(pipeline, "runRscript.sh")

configdict = TopmedPipeline.readConfig(configfile)
configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "data", "log"])

# analysis init
cluster.analysisInit(print_only=print_only)


job = "array_subset"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["subset_gds_file"] = configdict["subset_gds_file"] + ".chr .tmp"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile, version], array_range=chromosomes, mail=email, print_only=print_only)


job = "merge_gds"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["gds_file"] = configdict["subset_gds_file"] + ".chr .tmp"
config["merged_gds_file"] = configdict["subset_gds_file"]
config["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[jobid], email=email, print_only=print_only)


# remove temporary files (per-chr subset)
tmpid = cluster.submitJob(binary=True, job_name="rm_tmp_files", cmd="rm", args=[configdict["subset_gds_file"] + ".chr*.tmp"], holdid=[jobid], email=email, print_only=print_only)


# post analysis
job = "post_analysis"
jobpy = job + ".py"
pcmd=os.path.join(pipeline, jobpy)
argList = [pcmd, "-a", cluster.getAnalysisName(), "-l", cluster.getAnalysisLog(),
           "-s", cluster.getAnalysisStartSec()]
pdriver=os.path.join(pipeline, "run_python.sh")
cluster.submitJob(job_name=job, cmd=pdriver, args=argList, holdid=[jobid], print_only=print_only)
