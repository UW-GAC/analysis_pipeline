#! /usr/local/bin/python2.7

"""Convert VCF to GDS"""

import TopmedPipeline
import sys
import os
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
parser.add_argument("--clustertype", default="sge", 
                    help="type of compute cluster environment [default %(default)s]")
parser.add_argument("--clusterfile", default=None, 
                    help="file containing options to pass to the cluster (sge_request format)")
parser.add_argument("-c", "--chromosomes", default="1-22",
                    help="range of chromosomes [default %(default)s]")
parser.add_argument("-n", "--ncores", default="1-8",
                    help="number of cores to use; either a number (e.g, 1) or a range of numbers (e.g., 1-4) [default %(default)s]")
parser.add_argument("-e", "--email", default=None,
                    help="email address for job reporting")
parser.add_argument("--printOnly", action="store_true", default=False,
                    help="print qsub commands without submitting")
args = parser.parse_args()

configfile = args.configfile
chromosomes = args.chromosomes
clusterfile = args.clusterfile
clustertype = args.clustertype
ncores = args.ncores
email = args.email
printOnly = args.printOnly

opts = TopmedPipeline.getOptions(clusterfile)
cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type=clustertype, options=opts)

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))
driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
    
configdict = TopmedPipeline.readConfig(configfile)


job = "vcf2gds"

rscript = os.path.join(pipeline, "R", job + ".R")

# parsing bcf files relies on streaming bcftools output, so can't run in parallel
if os.path.splitext(configdict["vcf_file"])[1] == ".bcf":
    ncores = None

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile], array_range=chromosomes, request_cores=ncores, email=email, printOnly=printOnly)


job = "merge_gds"

rscript = os.path.join(pipeline, "R", job + ".R")

configdict["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(configdict, configfile)

holdid = [jobid["vcf2gds"]]

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], holdid=holdid, email=email, printOnly=printOnly)


job = "unique_variant_ids"

rscript = os.path.join(pipeline, "R", job + ".R")

holdid = [jobid["merge_gds"]]

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], holdid=holdid, email=email, printOnly=printOnly)
