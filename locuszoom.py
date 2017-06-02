#! /usr/local/bin/python2.7

"""LocusZoom"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
LocusZoom
"""
parser = ArgumentParser(description=description)
parser.add_argument("config_file", help="configuration file")
parser.add_argument("-c", "--chromosome", help="chromosome (1-24 or X,Y)")
parser.add_argument("--variant", help="variant ID")
parser.add_argument("--pop", help="population (AFR, AMR, ASN, EUR)")
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
args = parser.parse_args()

configfile = args.config_file
chromosome = args.chromosome
variant = args.variant
pop = args.pop
cluster_file = args.cluster_file
cluster_type = args.cluster_type
email = args.email
print_only = args.print_only
verbose = args.verbose

cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type, cluster_file, verbose)

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))
driver = os.path.join(pipeline, "runRscript.sh")

configdict = TopmedPipeline.readConfig(configfile)
configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "log", "plots"])


job = "locuszoom"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_prefix"] = configdict["plots_prefix"]
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

args = [rscript, configfile, "--chromosome " + chromosome, "--variant " + variant, "--pop " + pop]
jobid = cluster.submitJob(job_name=job, cmd=driver, args=args, email=email, print_only=print_only)


cluster.submitJob(job_name="cleanup", cmd=os.path.join(pipeline, "cleanup.sh"), holdid=jobid, print_only=print_only)
