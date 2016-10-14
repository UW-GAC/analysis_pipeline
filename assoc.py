#! /usr/local/bin/python2.7

"""Association tests"""

import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
Association tests
"""

parser = ArgumentParser(description=description)
parser.add_argument("assoctype", choices=["single", "window", "aggregate"],
                    help="type of association test")
parser.add_argument("configfile", help="configuration file")
parser.add_argument("-c", "--chromosomes", default="1-23",
                    help="range of chromosomes [default %(default)s]")
parser.add_argument("-p", "--pipeline", 
                    default="/projects/topmed/working_code/analysis_pipeline",
                    help="pipeline source directory")
parser.add_argument("-q", "--queue", default="olga.q", 
                    help="cluster queue name [default %(default)s]")
parser.add_argument("-e", "--email", default=None,
                    help="email address for job reporting")
parser.add_argument("--printOnly", action="store_true", default=False,
                    help="print qsub commands without submitting")
args = parser.parse_args()

assoctype = args.assoctype
configfile = args.configfile
chromosomes = args.chromosomes
pipeline = args.pipeline
queue = args.queue
email = args.email
printOnly = args.printOnly

sys.path.append(pipeline)
import TopmedPipeline
driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
    
configdict = TopmedPipeline.readConfig(configfile)

# check type of association test - single-variant unrelated is handled differently
single_unrel = assoctype == "single" and configdict["pcrelate_file"] == "NA"

if not single_unrel:
    job = "null_model"

    rscript = os.path.join(pipeline, "R", job + ".R")

    config = deepcopy(configdict)
    config["out_file"] = configdict["out_prefix"] + "_null_model.RData"
    configfile = configdict["out_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(config, configfile)

    jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], queue=queue, email=email, printOnly=printOnly)

    holdid = [jobid["null_model"]]

else:
    assoctype = "single_unrel"
    holdid = []


# for aggregate tests, generate variant list
if assoctype == "aggregate":
    job = "aggregate_list"
   
    rscript = os.path.join(pipeline, "R", job + ".R")

    config = deepcopy(configdict)
    config["out_file"] = configdict["out_prefix"] + "_" + job + "_chr .RData"
    configfile = configdict["out_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(config, configfile)

    jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], arrayRange=chromosomes, queue=queue, email=email, printOnly=printOnly)

    holdid.append(jobid["aggregate_list"].split(".")[0])


job = "assoc_" + assoctype

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["null_model_file"] = configdict["out_prefix"] + "_null_model.RData"
if assoctype == "aggregate":
    config["aggregate_variant_file"] = configdict["out_prefix"] + "_aggregate_list_chr .RData"
config["out_file"] = configdict["out_prefix"] + "_" + job + "_chr .RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], holdid=holdid, arrayRange=chromosomes, queue=queue, email=email, printOnly=printOnly)


prevjob = job
job = "assoc_plots"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["assoc_file"] = configdict["out_prefix"] + "_" + prevjob + "_chr .RData"
config["assoc_type"] = assoctype if not single_unrel else "single"
config["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
config["out_file_manh"] = configdict["out_prefix"] + "_manh.png"
config["out_file_qq"] = configdict["out_prefix"] + "_qq.png"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = [jobid["assoc_" + assoctype].split(".")[0]]

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], holdid=holdid, queue=queue, email=email, printOnly=printOnly)

