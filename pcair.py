#! /usr/local/bin/python2.7

"""Identity By Descent"""

import sys
import os
import subprocess
from optparse import OptionParser
from copy import deepcopy

usage = """%prog [options] config

PCA with the following steps:
1) Find unrelated sample set
2) Select SNPs with LD pruning using unrelated samples
3) PCA (using unrelated set, then project relatives)
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


job = "find_unrelated"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_related_file"] = configdict["out_prefix"] + "_related.RData"
config["out_unrelated_file"] = configdict["out_prefix"] + "_unrelated.RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], queue=qname, email=email, printOnly=printOnly)


job = "ld_pruning"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["sample_include_file"] = configdict["out_prefix"] + "_unrelated.RData"
config["out_file"] = configdict["out_prefix"] + "_pruned_variants_chr .RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = [jobid["find_unrelated"]]

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], holdid=holdid, arrayRange=chromosomes, queue=qname, email=email, printOnly=printOnly)


job = "combine_variants"

rscript = os.path.join(pipeline, "R", job + ".R")

config = dict()
chromRange = [int(x) for x in chromosomes.split("-")]
start = chromRange[0]
end = start if len(chromRange) == 1 else chromRange[1]
config["chromosomes"] = " ".join([str(x) for x in range(start, end + 1)])
config["in_file"] = configdict["out_prefix"] + "_pruned_variants_chr .RData"
config["out_file"] = configdict["out_prefix"] + "_pruned_variants.RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = [jobid["ld_pruning"].split(".")[0]]

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], holdid=holdid, queue=qname, email=email, printOnly=printOnly)


job = "pca_byrel"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["related_file"] = configdict["out_prefix"] + "_related.RData"
config["unrelated_file"] = configdict["out_prefix"] + "_unrelated.RData"
config["variant_include_file"] = configdict["out_prefix"] + "_pruned_variants.RData"
config["out_file"] = configdict["out_prefix"] + "_pcair.RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = [jobid["combine_variants"]]

jobid[job] = TopmedPipeline.submitJob(job, driver, [rscript, configfile], holdid=holdid, queue=qname, email=email, requestCores=multithread, printOnly=printOnly)

