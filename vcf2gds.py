#! /usr/local/bin/python2.7

"""Convert VCF to GDS"""

import sys
import os
import subprocess
from argparse import ArgumentParser
from copy import deepcopy

description = """
Convert VCF to GDS with the following steps:
1) Convert per-chromosome VCF files to GDS
2) Merge genotypes from per-chromosomes GDS files into a combined file
3) Assign unique variant id from merged file to per-chromosome GDS files
"""

parser = ArgumentParser(description=description)
parser.add_argument("configfile", help="configuration file")
parser.add_argument("-c", "--chromosomes", default="1-22",
                    help="range of chromosomes [default %(default)s]")
parser.add_argument("-p", "--pipeline", 
                    default="/projects/geneva/gcc-fs2/GCC_Code/analysis_pipeline",
                    help="pipeline source directory")
parser.add_argument("-q", "--queue", default="olga.q", 
                    help="cluster queue name [default %(default)s]")
parser.add_argument("-n", "--ncores", default="1-8",
                    help="number of cores to use; either a number (e.g, 1) or a range of numbers (e.g., 1-4) [default %(default)s]")
parser.add_argument("-e", "--email", default=None,
                    help="email address for job reporting")
parser.add_argument("--printOnly", action="store_true", default=False,
                    help="print qsub commands without submitting")
args = parser.parse_args()

configfile = args.configfile
chromosomes = args.chromosomes
pipeline = args.pipeline
queue = args.queue
ncores = args.ncores
email = args.email
printOnly = args.printOnly

sys.path.append(pipeline)
import TopmedPipeline
driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
    
configdict = TopmedPipeline.readConfig(configfile)


job = "vcf2gds"

rscript = os.path.join(pipeline, "R", job + ".R")

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], arrayRange=chromosomes, queue=queue, email=email, requestCores=ncores, printOnly=printOnly)


job = "merge_gds"

rscript = os.path.join(pipeline, "R", job + ".R")

chromRange = [int(x) for x in chromosomes.split("-")]
start = chromRange[0]
end = start if len(chromRange) == 1 else chromRange[1]
configdict["chromosomes"] = " ".join([str(x) for x in range(start, end + 1)])
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(configdict, configfile)

holdid = [jobid["vcf2gds"].split(".")[0]]

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], holdid=holdid, queue=queue, email=email, printOnly=printOnly)


job = "unique_variant_ids"

rscript = os.path.join(pipeline, "R", job + ".R")

holdid = [jobid["merge_gds"]]

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], holdid=holdid, queue=queue, email=email, printOnly=printOnly)
