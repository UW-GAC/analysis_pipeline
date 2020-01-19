#! /usr/bin/env python3

"""GRM"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy
from datetime import datetime, timedelta

description = """
Genetic Relationship Matrix (GRM)
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
                    help="print cluster commands without submitting")
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
submitPath = cluster.getSubmitPath()
driver = os.path.join(submitPath, "runRscript.sh")

configdict = TopmedPipeline.readConfig(configfile)
configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "data", "log"])

# analysis init
cluster.analysisInit(print_only=print_only)

job = "grm"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_file"] = configdict["data_prefix"] + "_grm_chr .gds"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile, version], request_cores=ncores, array_range=chromosomes, email=email, print_only=print_only)


job = "grm_combine"

rscript = os.path.join(pipeline, "R", job + ".R")

config = dict()
config["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
config["in_file"] = configdict["data_prefix"] + "_grm_chr .gds"
config["out_file"] = configdict["data_prefix"] + "_grm.gds"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[jobid], email=email, print_only=print_only)

# post analysis
bname = "post_analysis"
job = "grm" + "_" + bname
jobpy = bname + ".py"
pcmd=os.path.join(submitPath, jobpy)
argList = ["-a", cluster.getAnalysisName(), "-l", cluster.getAnalysisLog(),
           "-s", cluster.getAnalysisStartSec()]
cluster.submitJob(binary=True, job_name=job, cmd=pcmd, args=argList,
                  holdid=[jobid], print_only=print_only)
