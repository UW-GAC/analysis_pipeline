#! /usr/bin/env python3

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
parser.add_argument("--no_filter", action="store_true", default=False,
                    help="run VCF subsetting without filtering on >=1 alternate allele")
parser.add_argument("--md5", action="store_true", default=False,
                    help="generate md5 checksums")
parser.add_argument("--check_gds", action="store_true", default=False,
                    help="check resulting VCF against GDS")

args = parser.parse_args()

configfile = args.config_file
chromosomes = args.chromosomes
cluster_file = args.cluster_file
cluster_type = args.cluster_type
email = args.email
print_only = args.print_only
verbose = args.verbose
no_filter = args.no_filter
md5 = args.md5
check_gds = args.check_gds

version = "--version " + TopmedPipeline.__version__

cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type, cluster_file, verbose)

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))
submitPath = cluster.getSubmitPath()
driver = os.path.join(submitPath, "runRscript.sh")

configdict = TopmedPipeline.readConfig(configfile)
configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "log"])

# analysis init
cluster.analysisInit(print_only=print_only)

# submit job for each chromosome
job = "vcf_subset"
if no_filter:
    driver = os.path.join(submitPath, job + "_nofilter.sh")
else:
    driver = os.path.join(submitPath, job + ".sh")
samp_file = configdict["sample_file"]
in_file = configdict["vcf_file"].split(" ")
out_file = configdict["out_file"].split(" ")
args = [samp_file] + in_file + out_file
jobid = cluster.submitJob(job_name=job, cmd=driver, args=args, array_range=chromosomes, email=email, print_only=print_only)
holdid = [jobid]

# md5 checksums
if md5:
    job = "vcf_md5"
    driver = os.path.join(submitPath, job + ".sh")
    args = [os.path.dirname(configdict["out_file"])]
    
    md5id = cluster.submitJob(job_name=job, cmd=driver, args=args, holdid=[jobid], email=email, print_only=print_only)
    holdid.append(md5id)

# check file
if check_gds:
    job = "check_gds"
    driver = os.path.join(submitPath, "runRscript.sh")
    rscript = os.path.join(pipeline, "R", job + ".R")

    config = deepcopy(configdict)
    config["vcf_file"] = configdict["out_file"]
    configfile = configdict["config_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(config, configfile)

    args = ["-c", rscript, configfile]

    # add hold on anther array job, element-wise
    checkid = cluster.submitJob(binary=True, hold_array = jobid, job_name=job, cmd=driver, args=args, array_range=chromosomes, email=email, print_only=print_only)
    holdid.append(checkid)

# post analysis
bname = "post_analysis"
job = "vcf_subset" + "_" + bname
jobpy = bname + ".py"
pcmd=os.path.join(submitPath, jobpy)
argList = ["-a", cluster.getAnalysisName(), "-l", cluster.getAnalysisLog(),
           "-s", cluster.getAnalysisStartSec()]
cluster.submitJob(binary=True, job_name=job, cmd=pcmd, args=argList,
                  holdid=holdid, print_only=print_only)
