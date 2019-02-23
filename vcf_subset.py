#! /usr/bin/env python2.7

"""Generate subset of VCF files based on sample list"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """

"""

parser = ArgumentParser(description=description)
parser.add_argument("config_file", help="configuration file")
parser.add_argument("--cluster_type", default="UW_Cluster",
                    help="type of compute cluster environment [default %(default)s]")
parser.add_argument("--cluster_file", default=None,
                    help="json file containing options to pass to the cluster")
parser.add_argument("--verbose", action="store_true", default=False,
                    help="enable verbose output to help debug")
parser.add_argument("-c", "--chromosomes", default="1-23",
                    help="range of chromosomes [default %(default)s]")
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

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))
driver = os.path.join(pipeline, "runRscript.sh")

configdict = TopmedPipeline.readConfig(configfile)
configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "log"])

# analysis init
cluster.analysisInit(print_only=print_only)

# submit job for each chromosome
job = "vcf_subset"
driver = os.path.join(pipeline, job + ".sh")
samp_file = configdict["sample_file"]
in_file = configdict["vcf_file"].split(" ")
out_file = configdict["out_file"].split(" ")
args = [samp_file] + in_file + out_file
jobid = cluster.submitJob(job_name=job, cmd=driver, args=args, array_range=chromosomes, email=email, print_only=print_only)

# check file
job = "check_gds"
driver = os.path.join(pipeline, "runRscript.sh")
rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["vcf_file"] = configdict["out_file"]
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

args = ["-c", rscript, configfile]

# add special submit option to make one array job hold on anther array job, element-wise
job_cmd = cluster.clusterCfg["submit_cmd"]
subOpts = deepcopy(cluster.clusterCfg["submit_opts"])
subOpts["-hold_jid_ad"] = jobid
jobid = cluster.executeJobCmd(subOpts, job_cmd=job_cmd, job_name=job, cmd=driver, args=args, array_range=chromosomes, email=email, print_only=print_only)

# post analysis
job = "post_analysis"
jobpy = job + ".py"
pcmd=os.path.join(pipeline, jobpy)
argList = [pcmd, "-a", cluster.getAnalysisName(), "-l", cluster.getAnalysisLog(),
           "-s", cluster.getAnalysisStartSec()]
pdriver=os.path.join(pipeline, "run_python.sh")
cluster.submitJob(job_name=job, cmd=pdriver, args=argList,
                  holdid=[jobid], print_only=print_only)
