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
parser.add_argument("assoc_type", choices=["single", "window", "aggregate"],
                    help="type of association test")
parser.add_argument("config_file", help="configuration file")
parser.add_argument("-c", "--chromosomes", default="1-23",
                    help="range of chromosomes [default %(default)s]")
parser.add_argument("--segment_length", default=default_segment_length,
                    help="segment length in kb [default %(default)s]")
parser.add_argument("--n_segments", default=None,
                    help="number of segments for the entire genome (overrides segment_length)")
parser.add_argument("--cluster_type", default="UW_Cluster",
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

assoc_type = args.assoc_type
configfile = args.config_file
chromosomes = args.chromosomes
segment_length = args.segment_length
n_segments = args.n_segments
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
configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "data", "log", "plots", "report"])

# analysis init
cluster.analysisInit(print_only=print_only)
# hold is a list of submit IDs. A submit ID is a dict:
#     {jobname: [jobids]}
hold_null_agg = []

# copy parameter file for report
if "null_model_params" in configdict:
    paramfile = os.path.basename(configdict["config_prefix"]) + "_null_model.config.null_model.params"
    copyfile(configdict["null_model_params"], paramfile)


# for aggregate tests, generate variant list
if assoc_type == "aggregate":
    job = "aggregate_list"

    rscript = os.path.join(pipeline, "R", job + ".R")

    config = deepcopy(configdict)
    config["out_file"] = configdict["data_prefix"] + "_" + job + "_chr .RData"
    configfile = configdict["config_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(config, configfile)

    submitID = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile, version], array_range=chromosomes, email=email, print_only=print_only)
    hold_null_agg.append(submitID)


# define segments
if segment_length == default_segment_length and n_segments is None:
    build = configdict.setdefault("genome_build", "hg38")
    segment_file = os.path.join(submitPath, "segments_" + build + ".txt")
    print("Using default segment file for build " + build + " with segment_length " + default_segment_length + " kb")
else:
    job = "define_segments"

    rscript = os.path.join(submitPath, "R", job + ".R")

    config = deepcopy(configdict)
    config["out_file"] = configdict["config_prefix"] + "_segments.txt"
    configfile = configdict["config_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(config, configfile)

    segment_file = config["out_file"]

    # run and wait for results
    print("Defining segments...")
    log_file = job + "_" + strftime("%Y-%m-%d-%H-%M-%S", localtime()) + ".log"
    if n_segments is not None:
        cmd = [driver, rscript, configfile, "--n_segments " + n_segments]
    else:
        cmd = [driver, rscript, configfile, "--segment_length " + segment_length]
    cluster.runCmd(job_name=job, cmd=cmd, logfile=log_file)


# set up config for association test
config = deepcopy(configdict)
config["assoc_type"] = assoc_type

if assoc_type == "aggregate":
    config["aggregate_variant_file"] = configdict["data_prefix"] + "_aggregate_list_chr .RData"

assocScript = "assoc_" + assoc_type
config["out_prefix"] = configdict["data_prefix"] + "_" + assocScript
config["segment_file"] = segment_file
configfile = configdict["config_prefix"] + "_" + assocScript + ".config"
TopmedPipeline.writeConfig(config, configfile)


# get segments for each chromosome
chrom_list = TopmedPipeline.parseChromosomes(chromosomes).split(" ")
segment_list = TopmedPipeline.getChromSegments(segment_file, chrom_list)
segment_str = ["-".join([str(i) for i in s]) for s in segment_list]
segments = dict(list(zip(chrom_list, segment_str)))


# run association tests
hold_combine = []
for chromosome in chrom_list:
    job_assoc = assocScript + "_chr" + chromosome
    rscript = os.path.join(pipeline, "R", assocScript + ".R")
    args = ["-s", rscript, configfile, "--chromosome " + chromosome, version]
    # no email for jobs by segment
    submitID = cluster.submitJob(job_name=job_assoc, cmd=driver, args=args, holdid=hold_null_agg, array_range=segments[chromosome], print_only=print_only)

    combScript = "assoc_combine"
    job_comb = combScript + "_chr" + chromosome
    rscript = os.path.join(pipeline, "R", combScript + ".R")
    args = [rscript, configfile, "--chromosome " + chromosome, version]
    hold_assoc = [submitID]
    submitID = cluster.submitJob(job_name=job_comb, cmd=driver, args=args, holdid=hold_assoc, email=email, print_only=print_only)

    hold_combine.append(submitID)
# plots
job = "assoc_plots"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["assoc_file"] = configdict["data_prefix"] + "_" + assocScript + "_chr .RData"
config["assoc_type"] = assoc_type
config["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
config["out_file_manh"] = configdict["plots_prefix"] + "_manh.png"
config["out_file_qq"] = configdict["plots_prefix"] + "_qq.png"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

submitID = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=hold_combine, email=email, print_only=print_only)
hold_plots = [submitID]

# analysis report
job = "assoc_report"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_file"] = configdict["out_prefix"] + "_analysis_report"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

submitID = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=hold_plots, email=email, print_only=print_only)
hold_report = [submitID]

# post analysis
bname = "post_analysis"
job = "assoc" + "_" + assoc_type + "_" + bname
jobpy = bname + ".py"
pcmd=os.path.join(submitPath, jobpy)
argList = ["-a", cluster.getAnalysisName(), "-l", cluster.getAnalysisLog(),
           "-s", cluster.getAnalysisStartSec()]
cluster.submitJob(binary=True, job_name=job, cmd=pcmd, args=argList,
                  holdid=hold_report, print_only=print_only)
