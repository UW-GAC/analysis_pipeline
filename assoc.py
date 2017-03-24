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

parser = ArgumentParser(description=description)
parser.add_argument("assoc_type", choices=["single", "window", "aggregate"],
                    help="type of association test")
parser.add_argument("config_file", help="configuration file")
parser.add_argument("-c", "--chromosomes", default="1-23",
                    help="range of chromosomes [default %(default)s]")
parser.add_argument("--segment_length", default="10000",
                    help="segment length in kb [default %(default)s]")
parser.add_argument("--n_segments", default=None,
                    help="number of segments for the entire genome (overrides segment_length)")
parser.add_argument("--cluster_type", default="uw", 
                    help="type of compute cluster environment [default %(default)s]")
parser.add_argument("--cluster_file", default=None, 
                    help="file containing options to pass to the cluster (sge_request format)")
parser.add_argument("-e", "--email", default=None,
                    help="email address for job reporting")
parser.add_argument("--print_only", action="store_true", default=False,
                    help="print qsub commands without submitting")
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

opts = TopmedPipeline.getOptions(cluster_file)
cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type=cluster_type, options=opts)

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))
driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
    
configdict = TopmedPipeline.readConfig(configfile)

# check type of association test - single-variant unrelated is handled differently
no_pcrel = "pcrelate_file" not in configdict or configdict["pcrelate_file"] == "NA"
no_grm = "grm_file" not in configdict or configdict["grm_file"] == "NA"
single_unrel = assoc_type == "single" and no_pcrel and no_grm


# null model
if not single_unrel:
    job = "null_model"

    rscript = os.path.join(pipeline, "R", job + ".R")

    config = deepcopy(configdict)
    config["out_file"] = configdict["out_prefix"] + "_null_model.RData"
    configfile = configdict["out_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(config, configfile)

    opts = cluster.memoryOptions(job)

    jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], email=email, opts=opts, print_only=print_only)


    holdid = [jobid["null_model"]]
    assocScript = "assoc_" + assoc_type

else:
    holdid = []
    assocScript = "assoc_single_unrel"


# for aggregate tests, generate variant list
if assoc_type == "aggregate":
    job = "aggregate_list"
   
    rscript = os.path.join(pipeline, "R", job + ".R")

    config = deepcopy(configdict)
    config["out_file"] = configdict["out_prefix"] + "_" + job + "_chr .RData"
    configfile = configdict["out_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(config, configfile)

    opts = cluster.memoryOptions(job)

    jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile], array_range=chromosomes, email=email, opts=opts, print_only=print_only)

    holdid.append(jobid["aggregate_list"])

    
# define segments
job = "define_segments"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_file"] = configdict["out_prefix"] + "_segments.txt"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)
    
segment_file = config["out_file"]

# run and wait for results
print "Defining segments..."
log_file = open(job + "_" + strftime("%Y-%m-%d-%H-%M-%S", localtime()) + ".log", 'w')
if n_segments is not None:
    args = [driver, rscript, configfile, "--n_segments " + n_segments]
else:
    args = [driver, rscript, configfile, "--segment_length " + segment_length]
subprocess.check_call(args, stdout=log_file, stderr=log_file)
log_file.close()


# set up config for association test
config = deepcopy(configdict)
config["assoc_type"] = assoc_type
config["null_model_file"] = configdict["out_prefix"] + "_null_model.RData"
if assoc_type == "aggregate":
    config["aggregate_variant_file"] = configdict["out_prefix"] + "_aggregate_list_chr .RData"
config["out_prefix"] = configdict["out_prefix"] + "_" + assocScript
config["segment_file"] = segment_file
configfile = configdict["out_prefix"] + "_" + assocScript + ".config"
TopmedPipeline.writeConfig(config, configfile)


# get segments for each chromosome
chrom_list = TopmedPipeline.parseChromosomes(chromosomes).split(" ")
segment_list = TopmedPipeline.getChromSegments(segment_file, chrom_list)
segment_str = ["-".join([str(i) for i in s]) for s in segment_list]
segments = dict(zip(chrom_list, segment_str))


# run association tests
jobid_chrom = dict()
for chromosome in chrom_list:
    job_assoc = assocScript + "_chr" + chromosome
    rscript = os.path.join(pipeline, "R", assocScript + ".R")
    args = ["-s", rscript, configfile, "--chromosome " + chromosome]
    opts = cluster.memoryOptions("assoc")
    # no email for jobs by segment
    jobid[job_assoc] = cluster.submitJob(job_name=job_assoc, cmd=driver, args=args, holdid=holdid, array_range=segments[chromosome], opts=opts, print_only=print_only)

    combScript = "assoc_combine"
    job_comb = combScript + "_chr" + chromosome
    rscript = os.path.join(pipeline, "R", combScript + ".R")
    args = [rscript, configfile, "--chromosome " + chromosome]
    opts = cluster.memoryOptions("assoc_combine")
    holdid_comb = [jobid[job_assoc].split(".")[0]]
    jobid_chrom[job_comb] = cluster.submitJob(job_name=job_comb, cmd=driver, args=args, holdid=holdid_comb, email=email, opts=opts, print_only=print_only)
    
jobid.update(jobid_chrom)


# plots
job = "assoc_plots"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["assoc_file"] = configdict["out_prefix"] + "_" + assocScript + "_chr .RData"
config["assoc_type"] = assoc_type
config["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
config["out_file_manh"] = configdict["out_prefix"] + "_manh.png"
config["out_file_qq"] = configdict["out_prefix"] + "_qq.png"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = jobid_chrom.values()

opts = cluster.memoryOptions(job)

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], holdid=holdid, email=email, opts=opts, print_only=print_only)

