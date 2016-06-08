#! /usr/local/bin/python2.7

"""PC-Relate"""

import sys
import os
import subprocess
from optparse import OptionParser
from copy import deepcopy

usage = """%prog [options] config

PC-Relate
"""

parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/gcc-fs2/GCC_Code/analysis_pipeline",
                  help="pipeline source directory")
parser.add_option("-q", "--queue", dest="qname", default="olga.q", 
                  help="cluster queue name [default %default]")
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("--printOnly", dest="printOnly", action="store_true", default=False,
                  help="print qsub commands without submitting")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

configfile = args[0]
pipeline = options.pipeline
qname = options.qname
email = options.email
printOnly = options.printOnly

sys.path.append(pipeline)
import TopmedPipeline
driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
    
#configdict = TopmedPipeline.readConfig(configfile)


job = "pcrelate"

rscript = os.path.join(pipeline, "R", job + ".R")

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], queue=qname, email=email, printOnly=printOnly)
