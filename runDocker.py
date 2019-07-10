#! /usr/bin/env python
# using docker sdk, runs the topmed docker image to execute an analysis
# R function (via the R driver in the pipeline within docker container).
# This script is executed when a job is submitted to slurm and the compute
# nodes have docker installed
#
# slurm specific environment variables:
# job array (and tasks within array) related:
#    SLURM_ARRAY_TASK_COUNT  (number of tasks)
#    SLURM_ARRAY_TASK_ID     (index of task)
#    SLURM_ARRAY_TASK_MAX    (max number of tasks)
#    SLURM_ARRAY_TASK_MIN    (min index of task)
#    SLURM_ARRAY_TASK_STEP   (index step size)
#    SLURM_ARRAY_JOB_ID      (array's master job id)
# names
#    SLURM_CLUSTER_NAME     Name of the cluster on which the job is executing.
#    SLURM_JOB_NAME
#    SLURM_JOB_ID
#    SLURM_JOB_PARTITION
# job info
#    SLURM_CPUS_PER_TASK    Number of cpus/cores requested per task
#    SLURM_JOB_DEPENDENCY
#    SLURM_MEM_PER_NODE


import      sys
import      os
import      time
import      json
import      getpass
from        argparse import ArgumentParser
from        copy   import deepcopy

dockersdk = True
try:
    import  docker
except ImportError:
    global dockersdk
    dockersdk = False

# init globals
fileversion = '1.0'
msgErrPrefix = '>>> Error: '
msgInfoPrefix = '>>> Info: '
debugPrefix = '>>> Vebose: '
verbose = False

defDockerImage = "uwgac/topmed-roybranch"
defRunCmd = "/usr/local/analysis_pipeline/runRscript.sh"
defVolumes = "/projects:/projects"
defMemlimit = "8g"

slurmEnvDict = {
        "SLURM_ARRAY_TASK_COUNT": None,
        "SLURM_ARRAY_TASK_ID": None,
        "SLURM_ARRAY_TASK_MAX": None,
        "SLURM_ARRAY_TASK_MIN": None,
        "SLURM_ARRAY_TASK_STEP": None,
        "SLURM_ARRAY_JOB_ID": None,
        "SLURM_CLUSTER_NAME": None,
        "SLURM_JOB_NAME": None,
        "SLURM_JOB_ID": None,
        "SLURM_JOB_PARTITION": None,
        "SLURM_CPUS_PER_TASK": None,
        "SLURM_MEM_PER_NODE": None,
        "SLURM_JOB_DEPENDENCY": None
}

def pInfo(msg):
    tmsg=time.asctime()
    print(msgInfoPrefix+tmsg+": "+msg)

def pError(msg):
    tmsg=time.asctime()
    print(msgErrPrefix+tmsg+": "+msg)

def pDebug(msg):
    if verbose:
        tmsg=time.asctime()
        print(debugPrefix+tmsg+": "+msg)

def getSlurmEnv():
    slurmEnv = slurmEnvDict
    for key in slurmEnvDict.keys():
        slurmEnv[key] = os.getenv(key)
    return slurmEnv

def Summary(hdr):
    print(hdr)
    print("\tVersion: " + fileversion)
    print("\tDocker container run parameters:")
    print("\t\tDocker image: " + dockerimage)
    print("\t\tRun command: " + runcmd)
    print("\t\tRun cmd args: " + runargs)
    print("\t\tVolumes opts: " + str(volumes))
    print("\t\tEnvironment opts: " + str(environment))
    print("\t\tMemory limit: " + str(mem_limit))
    print("\t\tWorking dir: " + working_dir)
    print("\tSlurm environment:")
    print("\t\tCluster name: " + str(slurmEnv["SLURM_CLUSTER_NAME"]))
    print("\t\tPartition name: " + str(slurmEnv["SLURM_JOB_PARTITION"]))
    print("\t\tJob name: " + str(slurmEnv["SLURM_JOB_NAME"]))
    print("\t\tJob id: " + str(slurmEnv["SLURM_JOB_ID"]))
    print("\t\tJob dependency: " + str(slurmEnv["SLURM_JOB_DEPENDENCY"]))
    print("\t\tCPUs per task: " + str(slurmEnv["SLURM_CPUS_PER_TASK"]))
    print("\t\tArray job id: " + str(slurmEnv["SLURM_ARRAY_JOB_ID"]))
    print("\t\tArray task id: " + str(slurmEnv["SLURM_ARRAY_TASK_ID"]))
    print("\t\tArray task max: " + str(slurmEnv["SLURM_ARRAY_TASK_MAX"]))
    print("\t\tMemory per node: " + str(slurmEnv["SLURM_MEM_PER_NODE"]))
    print("\tVerbose: " + str(verbose))
    print("\tTest: " + str(test))
    print("\tDocker sdk installed? " + str(dockersdk))

# command line parser
parser = ArgumentParser(description = "Via docker sdk, create a docker container to execute an analysis job")
# docker container run: image and command
parser.add_argument( "--dockerimage", default = defDockerImage,
                     help = "Docker image [default: " + defDockerImage + "]" )
parser.add_argument( "--runcmd", default = defRunCmd,
                     help = "full path of pipeline run command [default: " + defRunCmd + "]" )
parser.add_argument( "--runargs", default = "", help = "Run cmd arguments" )
# docker container run: kwargs (run options)
parser.add_argument( "--volumes", default = defVolumes,
                     help = "Bind opts (/ld1:/rd1;/z1:/z2 etc)[default: " + defVolumes + "]" )
parser.add_argument( "--environment",
                     help = "Env opts (x=y;w=z etc)" )
parser.add_argument( "--mem_limit",
                     help = "Max memory (g) of container [default: unlimited]" )
parser.add_argument("--working_dir", help = "working directory [default: current working directory]")
# submit options explicitly passed (not available from env variable)
parser.add_argument( "--verbose", action="store_true", default = False,
                     help = "Verbose output [default: False]")
parser.add_argument( "--test", action="store_true", default = False,
                     help = "Test mode - don't execute [default: False]")

args = parser.parse_args()
dockerkwargs = {}
dockerimage = args.dockerimage
runcmd = args.runcmd
runargs = args.runargs
volumes = args.volumes
environment = args.environment
mem_limit = args.mem_limit
working_dir = args.working_dir
if working_dir == None:
    working_dir = os.getenv('PWD')
verbose = args.verbose
test = args.test
# get the slurm environment vars
slurmEnv = getSlurmEnv()

# kwargs (docker run options)
if mem_limit != None:
    dockerkwargs["mem_limit"] = str(mem_limit) + "g"
# working dir
dockerkwargs["working_dir"] = working_dir
# volumes
optdelim = ";"
volumes = volumes.rstrip(optdelim)
volumes = volumes.split(optdelim)
dockerkwargs["volumes"] = volumes
# environment
key = "environment"
if environment != None:
    elist = environment.rstrip(optdelim)
    elist = elist.split(optdelim)
    dockerkwargs[key] = elist
# environment - handle SLURM_CPUS_PER_TASK which changes to NSLOTS
if slurmEnv["SLURM_CPUS_PER_TASK"] != None:
    eslots = "NSLOTS=" + slurmEnv["SLURM_CPUS_PER_TASK"]
    if key in dockerkwargs.keys():
        dockerkwargs[key].append = eslots
    else:
        dockerkwargs[key] = [eslots]
# environment - set SGE_TASK_ID if an array job
if slurmEnv["SLURM_ARRAY_TASK_ID"] != None:
    etid = "SGE_TASK_ID=" + slurmEnv["SLURM_ARRAY_TASK_ID"]
    if key in dockerkwargs.keys():
        dockerkwargs[key].append = etid
    else:
        dockerkwargs[key] = [etid]
# misc
dockerkwargs["detach"] = True
dockerkwargs["remove"] = True

# create full docker run command
dockerFullCommand = runcmd + " " + runargs
if verbose:
    Summary("Summary of " + __file__)
pInfo("Docker container run parameters:")
print("\tImage: " + dockerimage)
print("\tCommand: " + dockerFullCommand)
print("\tKwargs: " + str(dockerkwargs))
pInfo("Slurm info")
print("\tCluster: " + str(slurmEnv["SLURM_CLUSTER_NAME"]))
print("\tPartition: " + str(slurmEnv["SLURM_JOB_PARTITION"]))
print("\tJob name: " + str(slurmEnv["SLURM_JOB_NAME"]))
print("\tJob id: " + str(slurmEnv["SLURM_JOB_ID"]))
print("\tJob dependency: " + str(slurmEnv["SLURM_JOB_DEPENDENCY"]))
print("\tCPUs per task: " + str(slurmEnv["SLURM_CPUS_PER_TASK"]))
if slurmEnv["SLURM_ARRAY_JOB_ID"] == None:
    print("\tJob is not an array job")
else:
    print("\tArray job id: " + str(slurmEnv["SLURM_ARRAY_JOB_ID"]))
    print("\tArray task index: " + str(slurmEnv["SLURM_ARRAY_TASK_ID"]))
    print("\tArray task max: " + str(slurmEnv["SLURM_ARRAY_TASK_MAX"]))


# if testing, just exit; else run docker container via sdk
if test:
    pInfo("Testing and not running docker")
else:
    if dockersdk:
        pDebug("Running docker container ...")
        if dockersdk:
            try:
                client = docker.from_env()
                dc = client.containers.run(dockerimage, command=dockerFullCommand, **dockerkwargs)
                og=dc.logs(stream=True)
                for line in og:
                   print(line.strip())
            except Exception as e:
                pError("Docker container exception: " + str(e))
                sys.exit(2)
        else:
            pError("Docker sdk not installed.")
            sys.exit(2)
        pInfo("Docker run completed.")
    else:
        pInfo("Docker sdk not installed; cannot run docker.")
        sys.exit(2)
