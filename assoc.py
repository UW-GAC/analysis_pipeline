#! /usr/local/bin/python2.7

"""Association tests"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
Association tests
"""

parser = ArgumentParser(description=description)
parser.add_argument("assocType", choices=["single", "window", "aggregate"],
                    help="type of association test")
parser.add_argument("configfile", help="configuration file")
parser.add_argument("-c", "--chromosomes", default="1-23",
                    help="range of chromosomes [default %(default)s]")
parser.add_argument("--clustertype", default="uw", 
                    help="type of compute cluster environment [default %(default)s]")
parser.add_argument("--clusterfile", default=None, 
                    help="file containing options to pass to the cluster (sge_request format)")
parser.add_argument("-e", "--email", default=None,
                    help="email address for job reporting")
parser.add_argument("--printOnly", action="store_true", default=False,
                    help="print qsub commands without submitting")
args = parser.parse_args()

assocType = args.assocType
configfile = args.configfile
chromosomes = args.chromosomes
clusterfile = args.clusterfile
clustertype = args.clustertype
email = args.email
printOnly = args.printOnly

opts = TopmedPipeline.getOptions(clusterfile)
cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type=clustertype, options=opts)

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))
driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
    
configdict = TopmedPipeline.readConfig(configfile)

# check type of association test - single-variant unrelated is handled differently
no_pcrel = "pcrelate_file" not in configdict or configdict["pcrelate_file"] == "NA"
no_grm = "grm_file" not in configdict or configdict["grm_file"] == "NA"
single_unrel = assocType == "single" and no_pcrel and no_grm

if not single_unrel:
    job = "null_model"

    rscript = os.path.join(pipeline, "R", job + ".R")

    config = deepcopy(configdict)
    config["out_file"] = configdict["out_prefix"] + "_null_model.RData"
    configfile = configdict["out_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(config, configfile)

    opts = cluster.memoryOptions(job)

    jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], email=email, opts=opts, printOnly=printOnly)


    holdid = [jobid["null_model"]]
    assocScript = "assoc_" + assocType

else:
    holdid = []
    assocScript = "assoc_single_unrel"


# for aggregate tests, generate variant list
if assocType == "aggregate":
    job = "aggregate_list"
   
    rscript = os.path.join(pipeline, "R", job + ".R")

    config = deepcopy(configdict)
    config["out_file"] = configdict["out_prefix"] + "_" + job + "_chr .RData"
    configfile = configdict["out_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(config, configfile)

    opts = cluster.memoryOptions(job)

    jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile], array_range=chromosomes, email=email, opts=opts, printOnly=printOnly)

    holdid.append(jobid["aggregate_list"])


segment_file = os.path.join(pipeline, "segments.txt")

config = deepcopy(configdict)
config["assoc_type"] = assocType
config["null_model_file"] = configdict["out_prefix"] + "_null_model.RData"
if assocType == "aggregate":
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

jobid_chrom = dict()
for chromosome in chrom_list:
    job_assoc = assocScript + "_chr" + chromosome
    rscript = os.path.join(pipeline, "R", assocScript + ".R")
    args = ["-s", rscript, configfile, "--chromosome " + chromosome]
    opts = cluster.memoryOptions("assoc")
    # no email for jobs by segment
    jobid[job_assoc] = cluster.submitJob(job_name=job_assoc, cmd=driver, args=args, holdid=holdid, array_range=segments[chromosome], opts=opts, printOnly=printOnly)

    combScript = "assoc_combine"
    job_comb = combScript + "_chr" + chromosome
    rscript = os.path.join(pipeline, "R", combScript + ".R")
    args = [rscript, configfile, "--chromosome " + chromosome]
    opts = cluster.memoryOptions("assoc_combine")
    holdid_comb = [jobid[job_assoc].split(".")[0]]
    jobid_chrom[job_comb] = cluster.submitJob(job_name=job_comb, cmd=driver, args=args, holdid=holdid_comb, email=email, opts=opts, printOnly=printOnly)
    
jobid.update(jobid_chrom)


job = "assoc_plots"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["assoc_file"] = configdict["out_prefix"] + "_" + assocScript + "_chr .RData"
config["assoc_type"] = assocType
config["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
config["out_file_manh"] = configdict["out_prefix"] + "_manh.png"
config["out_file_qq"] = configdict["out_prefix"] + "_qq.png"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = jobid_chrom.values()

opts = cluster.memoryOptions(job)

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], holdid=holdid, email=email, opts=opts, printOnly=printOnly)

