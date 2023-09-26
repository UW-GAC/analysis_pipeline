#! /usr/bin/env python3

"""Association tests"""

import TopmedPipeline
import sys
import os
import subprocess
from time import localtime, strftime
from argparse import ArgumentParser
from copy import deepcopy
from shutil import copyfile
from datetime import datetime, timedelta

description = """
Association tests
"""

default_segment_length = "10000"

parser = ArgumentParser(description=description)
parser.add_argument("config_file", help="configuration file")
parser.add_argument("-c", "--chromosomes", default="1-22",
                    help="range of chromosomes [default %(default)s]")
parser.add_argument("--cluster_type", default="GAC_Slurm_Cluster",
                    help="type of compute cluster environment [default %(default)s]")
parser.add_argument("--cluster_file", default=None,
                    help="json file containing cluster options")
parser.add_argument("-e", "--email", default=None,
                    help="email address for job reporting")
parser.add_argument("--print_only", action="store_true", default=False,
                    help="print qsub commands without submitting")
parser.add_argument("--verbose", action="store_true", default=False,
                    help="enable verbose output to help debug")
parser.add_argument("--version", action="version",
                    version="TopmedPipeline "+TopmedPipeline.__version__,
                    help="show the version number and exit")
args = parser.parse_args()

configfile = args.config_file
cluster_file = args.cluster_file
cluster_type = args.cluster_type
chromosomes = args.chromosomes
print_only = args.print_only
verbose = args.verbose
email = args.email

version = "--version " + TopmedPipeline.__version__

cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type, cluster_file, verbose)

pipeline = cluster.getPipelinePath()
submitPath = cluster.getSubmitPath()
driver = os.path.join(submitPath, "runRscript.sh")

configdict = TopmedPipeline.readConfig(configfile)
configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "data", "log", "report"])

# analysis init
cluster.analysisInit(print_only=print_only)

# calculate variant scores
job = "calc_variant_score"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_prefix"] = configdict["data_prefix"] + "_variant_score"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

# Submit by chromosome
submitID = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile, version], array_range=chromosomes, email=email, print_only=print_only)


# combine variant scores
job = "null_model_fast_scoreSE"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["variant_score_file"] = configdict["data_prefix"] + "_variant_score_chr .RData"
config["out_prefix"] = configdict["data_prefix"] + "_null_model_fast_scoreSE"
config["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

submitID = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[submitID], email=email, print_only=print_only)

# fast score approx report
# job = "null_model_fast_score_SE_report"
#
# rscript = os.path.join(pipeline, "R", job + ".R")
#
# config = deepcopy(configdict)
# config["out_prefix"] = configdict["out_prefix"] + "_null_model"
# configfile = configdict["config_prefix"] + "_" + job + ".config"
# TopmedPipeline.writeConfig(config, configfile)
#
# submitID = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[submitID], email=email, print_only=print_only)

# post analysis
bname = "post_analysis"
job = "null_model_fast_scoreSE" + "_" + bname
jobpy = bname + ".py"
pcmd=os.path.join(submitPath, jobpy)
argList = ["-a", cluster.getAnalysisName(), "-l", cluster.getAnalysisLog(),
           "-s", cluster.getAnalysisStartSec()]
cluster.submitJob(binary=True, job_name=job, cmd=pcmd, args=argList,
                  holdid=[submitID], print_only=print_only)
