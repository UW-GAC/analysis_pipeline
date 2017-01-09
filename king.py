#! /usr/local/bin/python2.7

"""Identity By Descent"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
Identity by Descent with the following steps:
1) Select SNPs with LD pruning
2) IBD calculations with KING-robust
"""

parser = ArgumentParser(description=description)
parser.add_argument("configfile", help="configuration file")
parser.add_argument("-c", "--chromosomes", default="1-22",
                    help="range of chromosomes [default %(default)s]")
parser.add_argument("--clustertype", default="uw", 
                    help="type of compute cluster environment [default %(default)s]")
parser.add_argument("--clusterfile", default=None, 
                    help="file containing options to pass to the cluster (sge_request format)")
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


job = "ld_pruning"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_file"] = configdict["out_prefix"] + "_pruned_variants_chr .RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile], array_range=chromosomes, email=email, printOnly=printOnly)


job = "combine_variants"

rscript = os.path.join(pipeline, "R", job + ".R")

config = dict()
config["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
config["in_file"] = configdict["out_prefix"] + "_pruned_variants_chr .RData"
config["out_file"] = configdict["out_prefix"] + "_pruned_variants.RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = [jobid["ld_pruning"]]

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile], holdid=holdid, email=email, printOnly=printOnly)


job = "ibd_king"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["variant_include_file"] = configdict["out_prefix"] + "_pruned_variants.RData"
config["out_file"] = configdict["out_prefix"] + "_ibd_king.RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = [jobid["combine_variants"]]

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], holdid=holdid, request_cores=ncores, email=email, printOnly=printOnly)


job = "kinship_plots"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["kinship_file"] = configdict["out_prefix"] + "_ibd_king.RData"
config["kinship_method"] = "king"
config["out_file_all"] = configdict["out_prefix"] + "_kinship_all.pdf"
config["out_file_cross"] = configdict["out_prefix"] + "_kinship_cross.pdf"
config["out_file_study"] = configdict["out_prefix"] + "_kinship_study.pdf"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = [jobid["ibd_king"]]

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], holdid=holdid, email=email, printOnly=printOnly)
