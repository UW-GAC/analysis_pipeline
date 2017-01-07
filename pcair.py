#! /usr/local/bin/python2.7

"""PC-AiR"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
PCA with the following steps:
1) Find unrelated sample set
2) Select SNPs with LD pruning using unrelated samples
3) PCA (using unrelated set, then project relatives)
"""

parser = ArgumentParser(description=description)
parser.add_argument("configfile", help="configuration file")
parser.add_argument("-c", "--chromosomes", default="1-22",
                    help="range of chromosomes [default %(default)s]")
parser.add_argument("--clustertype", default="sge", 
                    help="type of compute cluster environment [default %(default)s]")
parser.add_argument("--clusterfile", default=None, 
                    help="file containing options to pass to the cluster (sge_request format)")
parser.add_argument("-n", "--ncores", default="1-8",
                    help="number of cores to use; either a number (e.g, 1) or a range of numbers (e.g., 1-4) [default %(default)s]")
parser.add_argument("-e", "--email", default=None,
                    help="email address for job reporting")
parser.add_argument("--mem", action="store_true", default=False,
                    help="allocate memory for each job")
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


job = "find_unrelated"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_related_file"] = configdict["out_prefix"] + "_related.RData"
config["out_unrelated_file"] = configdict["out_prefix"] + "_unrelated.RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

opts = cluster.memoryOptions(job)

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], email=email, opts=opts, printOnly=printOnly)


job = "ld_pruning"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["sample_include_file"] = configdict["out_prefix"] + "_unrelated.RData"
config["out_file"] = configdict["out_prefix"] + "_pruned_variants_chr .RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = [jobid["find_unrelated"]]

opts = cluster.memoryOptions(job)

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile], holdid=holdid, array_range=chromosomes, email=email, opts=opts, printOnly=printOnly)


job = "combine_variants"

rscript = os.path.join(pipeline, "R", job + ".R")

config = dict()
config["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
config["in_file"] = configdict["out_prefix"] + "_pruned_variants_chr .RData"
config["out_file"] = configdict["out_prefix"] + "_pruned_variants.RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = [jobid["ld_pruning"]]

opts = cluster.memoryOptions(job)

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], holdid=holdid, email=email, opts=opts, printOnly=printOnly)


job = "pca_byrel"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["related_file"] = configdict["out_prefix"] + "_related.RData"
config["unrelated_file"] = configdict["out_prefix"] + "_unrelated.RData"
config["variant_include_file"] = configdict["out_prefix"] + "_pruned_variants.RData"
config["out_file"] = configdict["out_prefix"] + "_pcair.RData"
config["out_file_unrel"] = configdict["out_prefix"] + "_pcair_unrel.RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = [jobid["combine_variants"]]

opts = cluster.memoryOptions(job)

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], holdid=holdid, request_cores=ncores, email=email, opts=opts, printOnly=printOnly)


job = "pca_plots"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["pca_file"] = configdict["out_prefix"] + "_pcair.RData"
config["out_file_scree"] = configdict["out_prefix"] + "_pca_scree.pdf"
config["out_file_pc12"] = configdict["out_prefix"] + "_pca_pc12.pdf"
config["out_file_parcoord"] = configdict["out_prefix"] + "_pca_parcoord.pdf"
config["out_file_pairs"] = configdict["out_prefix"] + "_pca_pairs.png"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = [jobid["pca_byrel"]]

opts = cluster.memoryOptions(job)

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], holdid=holdid, email=email, opts=opts, printOnly=printOnly)


job = "pca_corr"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["pca_file"] = configdict["out_prefix"] + "_pcair_unrel.RData"
config["out_file"] = configdict["out_prefix"] + "_pcair_corr_chr .RData"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = [jobid["pca_byrel"]]

opts = cluster.memoryOptions(job)

# single core only
jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile], holdid=holdid, array_range=chromosomes, email=email, opts=opts, printOnly=printOnly)


job = "pca_corr_plots"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
config["corr_file"] = configdict["out_prefix"] + "_pcair_corr_chr .RData"
config["out_prefix"] = configdict["out_prefix"] + "_pcair_corr"
configfile = configdict["out_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = [jobid["pca_corr"]]

opts = cluster.memoryOptions(job)

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], holdid=holdid, email=email, opts=opts, printOnly=printOnly)
