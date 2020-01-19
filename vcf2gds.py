#! /usr/bin/env python3

"""Convert VCF to GDS"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
Convert VCF to GDS with the following steps:
1) Convert per-chromosome VCF files to GDS
2) Merge genotypes from per-chromosomes GDS files into a combined file
3) Assign unique variant id from merged file to per-chromosome GDS files
"""

parser = ArgumentParser(description=description)
parser.add_argument("config_file", help="configuration file")
parser.add_argument("--merge", action="store_true", default=False,
                    help="merge genotypes from per-chromosomes GDS files into a combined file")
parser.add_argument("--cluster_type", default="UW_Cluster",
                    help="type of compute cluster environment [default %(default)s]")
parser.add_argument("--cluster_file", default=None,
                    help="json file containing options to pass to the cluster")
parser.add_argument("--verbose", action="store_true", default=False,
                    help="enable verbose output to help debug")
parser.add_argument("-c", "--chromosomes", default="1-23",
                    help="range of chromosomes [default %(default)s]")
parser.add_argument("-n", "--ncores", default="1-8",
                    help="number of cores to use; either a number (e.g, 1) or a range of numbers (e.g., 1-4) [default %(default)s]")
parser.add_argument("-e", "--email", default=None,
                    help="email address for job reporting")
parser.add_argument("--print_only", action="store_true", default=False,
                    help="print qsub commands without submitting")
parser.add_argument("--version", action="version",
                    version="TopmedPipeline "+TopmedPipeline.__version__,
                    help="show the version number and exit")
args = parser.parse_args()

configfile = args.config_file
merge = args.merge
chromosomes = args.chromosomes
cluster_file = args.cluster_file
cluster_type = args.cluster_type
ncores = args.ncores
email = args.email
print_only = args.print_only
verbose = args.verbose

version = "--version " + TopmedPipeline.__version__

cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type, cluster_file, verbose)

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))
submitPath = cluster.getSubmitPath()
driver = os.path.join(submitPath, "runRscript.sh")

configdict = TopmedPipeline.readConfig(configfile)
configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "log"])

# analysis init
cluster.analysisInit(print_only=print_only)

job = "vcf2gds"

rscript = os.path.join(pipeline, "R", job + ".R")

# parsing bcf files relies on streaming bcftools output, so can't run in parallel
if os.path.splitext(configdict["vcf_file"])[1] == ".bcf":
    ncores = None

jobid = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile, version], array_range=chromosomes, request_cores=ncores, email=email, print_only=print_only)


job = "unique_variant_ids"

rscript = os.path.join(pipeline, "R", job + ".R")

jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[jobid], email=email, print_only=print_only)


job = "check_gds"

rscript = os.path.join(pipeline, "R", job + ".R")

jobid_chk = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile, version], holdid=[jobid], array_range=chromosomes, email=email, print_only=print_only)


# do we have more than one chromosome? if not, skip merge
if merge and len(TopmedPipeline.chromosomeRangeToList(chromosomes)) > 1:
    job = "merge_gds"

    rscript = os.path.join(pipeline, "R", job + ".R")

    configdict["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
    configfile = configdict["config_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(configdict, configfile)

    jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[jobid], email=email, print_only=print_only)


    job = "check_merged_gds"

    rscript = os.path.join(pipeline, "R", job + ".R")

    jobid = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile, version], holdid=[jobid], array_range=chromosomes, email=email, print_only=print_only)

# post analysis
bname = "post_analysis"
job = "vcf2gds" + "_" + bname
jobpy = bname + ".py"
pcmd=os.path.join(submitPath, jobpy)
argList = ["-a", cluster.getAnalysisName(), "-l", cluster.getAnalysisLog(),
           "-s", cluster.getAnalysisStartSec()]
cluster.submitJob(binary=True, job_name=job, cmd=pcmd, args=argList,
                  holdid=[jobid, jobid_chk], print_only=print_only)
