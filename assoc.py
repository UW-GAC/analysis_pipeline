#! /usr/local/bin/python2.7

"""Association tests"""

import TopmedPipeline
import sys
import os
import subprocess
from time import localtime, strftime
from argparse import ArgumentParser
from copy import deepcopy

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

cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type, cluster_file, verbose)

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))
driver = os.path.join(pipeline, "runRscript.sh")

configdict = TopmedPipeline.readConfig(configfile)
configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "data", "log", "plots", "report"])


# check type of association test - single-variant unrelated is handled differently
no_pcrel = "pcrelate_file" not in configdict or configdict["pcrelate_file"] == "NA"
no_grm = "grm_file" not in configdict or configdict["grm_file"] == "NA"
single_unrel = assoc_type == "single" and no_pcrel and no_grm

holdids = []

# null model
if not single_unrel:
    job = "null_model"

    rscript = os.path.join(pipeline, "R", job + ".R")

    config = deepcopy(configdict)
    config["out_file"] = configdict["data_prefix"] + "_null_model.RData"
    configfile = configdict["config_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(config, configfile)

    jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], email=email, print_only=print_only)

    holdids.append(jobid)
    assocScript = "assoc_" + assoc_type

else:
    assocScript = "assoc_single_unrel"


# for aggregate tests, generate variant list
if assoc_type == "aggregate":
    job = "aggregate_list"

    rscript = os.path.join(pipeline, "R", job + ".R")

    config = deepcopy(configdict)
    config["out_file"] = configdict["data_prefix"] + "_" + job + "_chr .RData"
    configfile = configdict["config_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(config, configfile)

    jobid = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile], array_range=chromosomes, email=email, print_only=print_only)
    holdids.append(jobid)


# define segments
if segment_length == default_segment_length and n_segments is None:
    build = configdict.setdefault("genome_build", "hg19")
    segment_file = os.path.join(pipeline, "segments_" + build + ".txt")
    print("Using default segment file for build " + build + " with segment_length " + default_segment_length + " kb")
else:
    job = "define_segments"

    rscript = os.path.join(pipeline, "R", job + ".R")

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
config["null_model_file"] = configdict["data_prefix"] + "_null_model.RData"
if assoc_type == "aggregate":
    config["aggregate_variant_file"] = configdict["data_prefix"] + "_aggregate_list_chr .RData"
config["out_prefix"] = configdict["data_prefix"] + "_" + assocScript
config["segment_file"] = segment_file
configfile = configdict["config_prefix"] + "_" + assocScript + ".config"
TopmedPipeline.writeConfig(config, configfile)


# get segments for each chromosome
chrom_list = TopmedPipeline.parseChromosomes(chromosomes).split(" ")
segment_list = TopmedPipeline.getChromSegments(segment_file, chrom_list)
segment_str = ["-".join([str(i) for i in s]) for s in segment_list]
segments = dict(zip(chrom_list, segment_str))


# run association tests
holdids_combine = []
for chromosome in chrom_list:
    job_assoc = assocScript + "_chr" + chromosome
    rscript = os.path.join(pipeline, "R", assocScript + ".R")
    args = ["-s", rscript, configfile, "--chromosome " + chromosome]
    # no email for jobs by segment
    jobid = cluster.submitJob(job_name=job_assoc, cmd=driver, args=args, holdid=holdids, array_range=segments[chromosome], print_only=print_only)

    combScript = "assoc_combine"
    job_comb = combScript + "_chr" + chromosome
    rscript = os.path.join(pipeline, "R", combScript + ".R")
    args = [rscript, configfile, "--chromosome " + chromosome]
    jobid = cluster.submitJob(job_name=job_comb, cmd=driver, args=args, holdid=jobid, email=email, print_only=print_only)
    holdids_combine.append(jobid)


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

jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], holdid=holdids_combine, email=email, print_only=print_only)


# analysis report
job = "assoc_report"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_file"] = configdict["out_prefix"] + "_analysis_report"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], holdid=jobid, email=email, print_only=print_only)


cluster.submitJob(job_name="cleanup", cmd=os.path.join(pipeline, "cleanup.sh"), holdid=jobid, print_only=print_only)
