#! /usr/local/bin/python2.7

"""Association tests (single-variant)"""

import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
Association tests (single-variant)
"""

parser = ArgumentParser(description=description)
parser.add_argument("configfile", help="configuration file")
parser.add_argument("assoctype", choices=["single", "window"],
                    help="type of association test")
parser.add_argument("-c", "--chromosomes", default="1-23",
                    help="range of chromosomes [default %(default)s]")
parser.add_argument("-p", "--pipeline", 
                    default="/projects/geneva/gcc-fs2/GCC_Code/analysis_pipeline",
                    help="pipeline source directory")
parser.add_argument("-q", "--queue", default="olga.q", 
                    help="cluster queue name [default %(default)s]")
parser.add_argument("-e", "--email", default=None,
                    help="email address for job reporting")
parser.add_argument("--printOnly", action="store_true", default=False,
                    help="print qsub commands without submitting")
args = parser.parse_args()

configfile = args.configfile
assoctype = args.assoctype
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


job = "null_model"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_file"] = configdict["out_prefix"] + "_null_model.RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], queue=queue, email=email, printOnly=printOnly)


job = "assoc_" + assoctype

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["null_model_file"] = configdict["out_prefix"] + "_null_model.RData"
config["out_file"] = configdict["out_prefix"] + "_" + job + "_chr .RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = [jobid["null_model"]]

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], holdid=holdid, arrayRange=chromosomes, queue=queue, email=email, printOnly=printOnly)
