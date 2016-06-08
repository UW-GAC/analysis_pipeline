#! /usr/local/bin/python2.7

"""Convert VCF to GDS"""

import sys
import os
import subprocess
from optparse import OptionParser
from copy import deepcopy

usage = """%prog [options] config

Convert VCF to GDS with the following steps:
1) Convert per-chromosome VCF files to GDS
2) Merge genotypes from per-chromosomes GDS files into a combined file
3) Assign unique variant id from merged file to per-chromosome GDS files
"""

parser = OptionParser(usage=usage)
parser.add_option("-c", "--chromosomes", dest="chromosomes", default="1-22",
                  help="range of chromosomes [default %default]")
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/gcc-fs2/GCC_Code/analysis_pipeline",
                  help="pipeline source directory")
parser.add_option("-q", "--queue", dest="qname", default="olga.q", 
                  help="cluster queue name [default %default]")
parser.add_option("-m", "--multithread", dest="multithread", default="1-8",
                  help="number of cores to use; either a number (e.g, 1) or a range of numbers (e.g., 1-4) [default %default]")
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("--printOnly", dest="printOnly", action="store_true", default=False,
                  help="print qsub commands without submitting")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

configfile = args[0]
chromosomes = options.chromosomes
pipeline = options.pipeline
qname = options.qname
multithread = options.multithread
email = options.email
printOnly = options.printOnly

sys.path.append(pipeline)
import TopmedPipeline
driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
    
configdict = TopmedPipeline.readConfig(configfile)


job = "vcf2gds"

rscript = os.path.join(pipeline, "R", job + ".R")

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], arrayRange=chromosomes, queue=qname, email=email, requestCores=multithread, printOnly=printOnly)


job = "merge_gds"

rscript = os.path.join(pipeline, "R", job + ".R")

chromRange = [int(x) for x in chromosomes.split("-")]
start = chromRange[0]
end = start if len(chromRange) == 1 else chromRange[1]
configdict["chromosomes"] = " ".join([str(x) for x in range(start, end + 1)])
configfile = job + ".config"
TopmedPipeline.writeConfig(configdict, configfile)

holdid = [jobid["vcf2gds"].split(".")[0]]

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], holdid=holdid, queue=qname, email=email, printOnly=printOnly)


job = "unique_variant_ids"

rscript = os.path.join(pipeline, "R", job + ".R")

holdid = [jobid["merge_gds"]]

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], holdid=holdid, queue=qname, email=email, printOnly=printOnly)
