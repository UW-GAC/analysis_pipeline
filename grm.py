#! /usr/local/bin/python2.7

"""GRM"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
PC-Relate
"""
parser = ArgumentParser(description=description)
parser.add_argument("config_file", help="configuration file")
parser.add_argument("--cluster_type", default="uw", 
                    help="type of compute cluster environment [default %(default)s]")
parser.add_argument("--cluster_file", default=None, 
                    help="file containing options to pass to the cluster (sge_request format)")
parser.add_argument("-n", "--ncores", default="1-8",
                    help="number of cores to use; either a number (e.g, 1) or a range of numbers (e.g., 1-4) [default %(default)s]")
parser.add_argument("-e", "--email", default=None,
                    help="email address for job reporting")
parser.add_argument("--print_only", action="store_true", default=False,
                    help="print cluster commands without submitting")
args = parser.parse_args()

configfile = args.config_file
cluster_file = args.cluster_file
cluster_type = args.cluster_type
ncores = args.ncores
email = args.email
print_only = args.print_only

opts = TopmedPipeline.getOptions(cluster_file)
cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type=cluster_type, options=opts)

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))
driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
    
configdict = TopmedPipeline.readConfig(configfile)


job = "grm"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_file"] = configdict["out_prefix"] + "_grm.RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], request_cores=ncores, email=email, print_only=print_only)
