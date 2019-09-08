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
import      errno
import      os
import      time
import      json
import      getpass
import      subprocess
from        argparse import ArgumentParser
from        copy   import deepcopy

try:
    import  docker
except ImportError:
    dockersdk = False
else:
    dockersdk = True
# init globals
fileversion = '1.0'
msgErrPrefix = '>>> Error: '
msgInfoPrefix = '>>> Info: '
debugPrefix = '>>> Vebose: '
verbose = False

defDockerImage = "uwgac/topmed-roybranch:latest"
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
def GetSlurmJobID():
    # if array job, use the array job id; else job id
    # check if array_job
    sj = {"arrayjob":None, "jobid": None}
    if slurmEnv["SLURM_ARRAY_JOB_ID"] != None:
        sj["arrayjob"] = True
        sj["jobid"] = slurmEnv["SLURM_ARRAY_JOB_ID"]
    else:
        sj["arrayjob"] = False
        if slurmEnv["SLURM_JOB_ID"] != None:
            sj["jobid"] =  slurmEnv["SLURM_JOB_ID"]
        else:
            sj["jobid"] = "NOJOBID"
    return sj

def CreateLogFileName():
    # general name is <jobname>_<array_job_id>_<array_index>.log or
    #                 <jobname>_<job_id>.log
    ext = ".log"
    # get slurm job id
    slurmjid = GetSlurmJobID()
    jid = slurmjid["jobid"]
    # jobname
    if slurmEnv["SLURM_JOB_NAME"] != None:
        jn = slurmEnv["SLURM_JOB_NAME"]
    else:
        jn = "NOJOBNAME"
    # check if array_job
    if slurmjid["arrayjob"]:
        if slurmEnv["SLURM_ARRAY_TASK_ID"] != None:
            tid = slurmEnv["SLURM_ARRAY_TASK_ID"]
        else:
            tid = "NOTASKID"
        lfn = jn + "_" + jid + "_" + tid + ext
    else:
        lfn = jn + "_" + jid + ext
    return lfn


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
    print("\t\tDocker run command:\n\t\t" + dockerFullCommand)
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
    if log:
        print("\tLog file: " + logfile)
    print("\tVerbose: " + str(verbose))
    print("\tTest: " + str(test))
    print("\tDocker sdk installed? " + str(dockersdk))
# start the timer
startTime = time.time()
# command line parser
parser = ArgumentParser(description = "Via python Popen, run docker to execute an analysis job")
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
parser.add_argument("--machine", help = "machine type of partition")
parser.add_argument("--cost", help="hourly cost of machine type")
parser.add_argument("--log", action="store_true", default = False, help = "create a docker log file" )
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
machine = args.machine
cost = args.cost
verbose = args.verbose
test = args.test
log = args.log
# get the slurm environment vars
slurmEnv = getSlurmEnv()
# log
if log:
    logfile = CreateLogFileName()

# docker run options
dockeropts = ""
# standard opt (--rm)
dockeropts += "-a stdout -a stderr --rm "
if mem_limit != None:
    dockeropts += "-m " + str(mem_limit) + "g "
# working dir
dockeropts += "-w " + working_dir + " "
# volumes
optdelim = ";"
volumes = volumes.rstrip(optdelim)
volumes = volumes.split(optdelim)
for v in volumes:
    dockeropts += "-v " + v + " "
# environment
if environment != None:
    elist = environment.rstrip(optdelim)
    elist = elist.split(optdelim)
    for e in elist:
        dockeropts += "-e " + e + " "
# environment - set JOB_ID for runRscript.sh (in case of error)
slurmjid = GetSlurmJobID()
rjid = "JOB_ID=" + slurmjid["jobid"]
dockeropts += "-e " + rjid + " "
# environment - set SGE_TASK_ID if an array job
if slurmjid["arrayjob"]:
    sgeid = "SGE_TASK_ID=" + slurmEnv["SLURM_ARRAY_TASK_ID"]
    dockeropts += "-e " + sgeid + " "
# environment - handle SLURM_CPUS_PER_TASK which changes to NSLOTS
if slurmEnv["SLURM_CPUS_PER_TASK"] != None:
    slots = "NSLOTS=" + slurmEnv["SLURM_CPUS_PER_TASK"]
    dockeropts += "-e " + slots + " "
# create full docker run command
dockerFullCommand = "docker run " + dockeropts + dockerimage + " " + runcmd + " " + runargs
if log:
    dockerFullCommand +=  " > " + logfile + " 2>&1"
if verbose:
    Summary("Summary of " + __file__)
pInfo("Docker container run parameters:")
print("\tImage: " + dockerimage)
print("\tRun opts:")
print("\tCommand: " + dockerFullCommand)
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
    status = 0
else:
    if dockersdk:
        if machine != None:
            pInfo("Machine type: " + machine + " ( $" + str(cost) + "/hr )")
        try:
            pDebug("Getting docker client ...")
            dclient = docker.from_env()
            pDebug("Pulling latest ...")
            dimage = dclient.images.pull(dockerimage)
            cmd = dockerFullCommand
            cmdl = cmd.split()
            pDebug("cmdl = " + str(cmdl))
            pInfo("Executing docker via popen:\n\t" + cmd + " ...")
            if log:
                pInfo("Sending stdout/stderr of docker run to: " + logfile)
            process = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr, shell=True)
            status = process.wait()
            if status:
                smsg = os.strerror(status)
                if status == -9:
                    smsg = "possibly killed externally via kill -9"
                if status == 137:
                    smsg = "possibly killed internally via memory limit"
                pError("Executing docker failed. Error: " + str(status) + " - " + smsg)
            else:
                pInfo( "Executing docker completed without errors.")
        except Exception as e:
            pError("Docker container exception: " + str(e))
            status = 2
        eTime = time.time() - startTime
        eTimeHr = eTime/60./60.
        pInfo("Elapsed time (hr): " + str(eTimeHr))
        if status == 0:
            if cost != None:
                totalCost = eTimeHr*float(cost)
                pInfo("Estimated cost= " + "$" + str(totalCost))
            pInfo("Docker run completed successfully.")
        else:
            if cost != None:
                totalCost = eTimeHr*float(cost)
                pInfo("Estimated cost up to error= " + "$" + str(totalCost))
    else:
        pError("Docker sdk not installed; cannot run docker.")
        status = 2
# exit
    sys.exit(status)
