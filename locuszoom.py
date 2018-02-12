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

cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type, cluster_file, verbose)

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))
driver = os.path.join(pipeline, "runRscript.sh")

configdict = TopmedPipeline.readConfig(configfile)
configdict = TopmedPipeline.directorySetup(configdict, subdirs=["log", "plots"])


job = "locuszoom"

rscript = os.path.join(pipeline, "R", job + ".R")

# find number of jobs to submit by counting lines in file
n = TopmedPipeline.countLines(configdict["locus_file"])
range = "1-" + str(n-1)

args = ["-s", rscript, configfile]
jobid = cluster.submitJob(job_name=job, cmd=driver, args=args, array_range=range, email=email, print_only=print_only)


cluster.submitJob(job_name="cleanup", cmd=os.path.join(pipeline, "cleanup.sh"), holdid=jobid, print_only=print_only)
