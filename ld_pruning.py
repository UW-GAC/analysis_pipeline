#! /usr/bin/env python3

"""LD Pruning"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
LD pruning with the following steps:
1) Select variants with LD pruning
2) Combine selected variants from all chromosomes
3) Create GDS file with only pruned variants
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
submitPath = cluster.getSubmitPath()
driver = os.path.join(submitPath, "runRscript.sh")

configdict = TopmedPipeline.readConfig(configfile)
configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "data", "log", "plots"])

# analysis init
cluster.analysisInit(print_only=print_only)


job = "ld_pruning"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_file"] = configdict["data_prefix"] + "_pruned_variants_chr .RData"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile, version], array_range=chromosomes, email=email, print_only=print_only)


job = "subset_gds"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["subset_gds_file"] = configdict["subset_gds_file"] + ".chr .tmp"
config["variant_include_file"] = configdict["data_prefix"] + "_pruned_variants_chr .RData"
config["sample_include_file"] = "NA"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile, version], array_range=chromosomes, holdid=[jobid], email=email, print_only=print_only)


job = "merge_gds"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["gds_file"] = configdict["subset_gds_file"] + ".chr .tmp"
config["merged_gds_file"] = configdict["subset_gds_file"]
config["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[jobid], email=email, print_only=print_only)


job = "check_merged_gds"

rscript = os.path.join(pipeline, "R", job + ".R")

jobid = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile, version], holdid=[jobid], array_range=chromosomes, email=email, print_only=print_only)


# remove temporary files (per-chr subset)
tmpid = cluster.submitJob(binary=True, job_name="rm_tmp_files", cmd="rm", args=[configdict["subset_gds_file"] + ".chr*.tmp"], holdid=[jobid], email=email, print_only=print_only)

# post analysis
bname = "post_analysis"
job = "ld_pruning" + "_" + bname
jobpy = bname + ".py"
pcmd=os.path.join(submitPath, jobpy)
argList = ["-a", cluster.getAnalysisName(), "-l", cluster.getAnalysisLog(),
           "-s", cluster.getAnalysisStartSec()]
cluster.submitJob(binary=True, job_name=job, cmd=pcmd, args=argList,
                  holdid=[jobid], print_only=print_only)
