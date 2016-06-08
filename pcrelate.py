#! /usr/local/bin/python2.7

"""PC-Relate"""

import sys
import os
import subprocess
from argparse import ArgumentParser
from copy import deepcopy

description = """
PC-Relate
"""
parser = ArgumentParser(description=description)
parser.add_argument("configfile", help="configuration file")
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
pipeline = args.pipeline
queue = args.queue
email = args.email
printOnly = args.printOnly

sys.path.append(pipeline)
import TopmedPipeline
driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
    
#configdict = TopmedPipeline.readConfig(configfile)


job = "pcrelate"

rscript = os.path.join(pipeline, "R", job + ".R")

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], queue=queue, email=email, printOnly=printOnly)
