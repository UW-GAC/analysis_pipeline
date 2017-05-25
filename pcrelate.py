#! /usr/local/bin/python2.7

"""PC-Relate"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
PC-Relate
"""
parser = ArgumentParser(description=description)
parser.add_argument("config_file", help="configuration file")
parser.add_argument("--cluster_type", default="uw", 
                    help="type of compute cluster environment [default %(default)s]")
parser.add_argument("--cluster_file", default=None, 
                    help="file containing options to pass to the cluster (sge_request format)")
parser.add_argument("-e", "--email", default=None,
                    help="email address for job reporting")
parser.add_argument("--print_only", action="store_true", default=False,
                    help="print cluster commands without submitting")
args = parser.parse_args()

configfile = args.config_file
cluster_file = args.cluster_file
cluster_type = args.cluster_type
email = args.email
print_only = args.print_only

opts = TopmedPipeline.getOptions(cluster_file)
cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type=cluster_type, options=opts)

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))
driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
    
configdict = TopmedPipeline.readConfig(configfile)
configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "data", "log", "plots"])


job = "pcrelate"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_prefix"] = configdict["data_prefix"]
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], email=email, print_only=print_only)


job = "kinship_plots"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["kinship_file"] = configdict["data_prefix"] + "_pcrelate.gds"
config["kinship_method"] = "pcrelate"
config["out_file_all"] = configdict["plots_prefix"] + "_kinship_all.pdf"
config["out_file_cross"] = configdict["plots_prefix"] + "_kinship_cross.pdf"
config["out_file_study"] = configdict["plots_prefix"] + "_kinship_study.pdf"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

holdid = [jobid["pcrelate"]]

jobid[job] = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile], holdid=holdid, email=email, print_only=print_only)


cluster.submitJob(job_name="cleanup", cmd=os.path.join(pipeline, "cleanup.sh"), holdid=[jobid["kinship_plots"]], print_only=print_only)
