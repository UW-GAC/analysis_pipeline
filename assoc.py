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
parser.add_argument("-q", "--queue", default="olga.q", 
                    help="cluster queue name [default %(default)s]")
parser.add_argument("-e", "--email", default=None,
                    help="email address for job reporting")
parser.add_argument("--printOnly", action="store_true", default=False,
                    help="print qsub commands without submitting")
args = parser.parse_args()

assocType = args.assocType
configfile = args.configfile
chromosomes = args.chromosomes
queue = args.queue
email = args.email
printOnly = args.printOnly

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))
driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
    
configdict = TopmedPipeline.readConfig(configfile)

qsubOpts = ""

# check type of association test - single-variant unrelated is handled differently
single_unrel = assocType == "single" and configdict["pcrelate_file"] == "NA"

if not single_unrel:
    job = "null_model"

    rscript = os.path.join(pipeline, "R", job + ".R")

    config = deepcopy(configdict)
    config["out_file"] = configdict["out_prefix"] + "_null_model.RData"
    configfile = configdict["out_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(config, configfile)

    #qsubOpts = "-l h_vmem=1.2G"
    jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], queue=queue, email=email, qsubOptions=qsubOpts, printOnly=printOnly)

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

    #qsubOpts = "-l h_vmem=3G"
    jobid[job] = TopmedPipeline.submitJob(job, driver, ["-c", rscript, configfile], arrayRange=chromosomes, queue=queue, email=email, qsubOptions=qsubOpts, printOnly=printOnly)

    holdid.append(jobid["aggregate_list"].split(".")[0])


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

#jobid[assocScript] = TopmedPipeline.submitJob(assocScript, driver, ["-c", rscript, configfile], holdid=holdid, arrayRange=chromosomes, queue=queue, email=email, printOnly=printOnly)

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
    # no email for jobs by segment
    #qsubOpts = "-l h_vmem=3G"
    jobid[job_assoc] = TopmedPipeline.submitJob(job_assoc, driver, args, holdid=holdid, arrayRange=segments[chromosome], queue=queue, qsubOptions=qsubOpts, printOnly=printOnly)

    combScript = "assoc_combine"
    job_comb = combScript + "_chr" + chromosome
    rscript = os.path.join(pipeline, "R", combScript + ".R")
    args = [rscript, configfile, "--chromosome " + chromosome]
    holdid_comb = [jobid[job_assoc].split(".")[0]]
    jobid_chrom[job_comb] = TopmedPipeline.submitJob(job_comb, driver, args, holdid=holdid_comb, queue=queue, email=email, qsubOptions=qsubOpts, printOnly=printOnly)
    
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

#holdid = [jobid[assocScript].split(".")[0]]
holdid = [c.split(".")[0] for c in jobid_chrom.values()]

#qsubOpts = "-l h_vmem=3.2G"
jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], holdid=holdid, queue=queue, email=email, qsubOptions=qsubOpts, printOnly=printOnly)

