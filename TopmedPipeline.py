#! /usr/bin/env python3
from __future__ import division
"""Utility functions for TOPMed pipeline"""

__version__ = "2.8.1"

import os
import sys
import csv
import subprocess
from copy import deepcopy
import getpass
import time
import json
import math
import collections
import datetime
import awsbatch
import port_popen

try:
    import boto3
except ImportError:
    batchSupport = False
else:
    batchSupport = True



def readConfig(file):
    """Read a pipeline config file.

    Usage:
    config = readConfig(file)

    Arguments:
    file - name of config file to read

    Returns:
    dictionary with config values
    """

    config = {}
    f = open(file, 'r')
    reader = csv.reader(f, delimiter=' ', quotechar='"', skipinitialspace=True)
    for line in reader:
        if len(line) == 0:
            continue

        if line[0][0] == "#":
            continue

        if len(line) == 1:
            sys.exit("Error reading config file " + file + ":\nNo value for parameter '" + line[0] + "' in line " + str(reader.line_num))

        if len(line) > 2:
            if line[2] == '':
                line = line[0:2]
            else:
                sys.exit("Error reading config file " + file + ":\nToo many parameters in line " + str(reader.line_num))

        (key, value) = line
        config[key] = value

    f.close()
    return config


def writeConfig(config, file):
    """Write a pipeline config file.

    Usage:
    writeConfig(config, file)

    Arguments:
    config - dict with config parameters
    file - name of config file to write
    """

    f = open(file, 'w')
    writer = csv.writer(f, delimiter=' ', quotechar='"')
    for key, value in config.items():
        writer.writerow([key, value])
    f.close()


def getFirstColumn(file, skipHeader=True):
    """Read a file and return the first column

    Usage:
    x = getFirstColumn(file)

    Arguments:
    file - name of file to read

    Returns:
    list with values in the first column (minus the header)
    """
    f = open(file, 'r')
    reader = csv.reader(f, delimiter="\t")
    if skipHeader:
        dummy = next(reader)
    x = [line[0] for line in reader]
    f.close()

    return x


def which(x, y):
    """Returns indices of x that equal y (1-based)
    """
    return [ i+1 for i, j in enumerate(x) if j == y ]


def getChromSegments(map_file, chromosome):
    """Read a pipeline segments file.

    Usage:
    segments = getChromSegments(map_file, chromosome)

    Arguments:
    file - name of segments file to read (expect first column is chromosome)
    chromosome - character value for chromosome

    Returns:
    list with beginning and ending segment indices for each chromosome
    """
    chrom_segments = getFirstColumn(map_file)

    # get indices of segments matching this chromosome
    segments = [ (min(x), max(x)) for x in [ which(chrom_segments, c) for c in chromosome ] ]

    return segments


def chromosomeRangeToList(chromosomes):
    chromRange = [int(x) for x in chromosomes.split("-")]
    start = chromRange[0]
    end = start if len(chromRange) == 1 else chromRange[1]
    return list(range(start, end + 1))

def parseChromosomes(chromosomes):
    chromString = " ".join([str(x) for x in chromosomeRangeToList(chromosomes)])
    chromString = chromString.replace("23", "X")
    chromString = chromString.replace("24", "Y")
    return chromString


def dictToString(d):
    """Construct a string from a dictionary"""
    s = ' '.join([k + ' ' + v for k, v in d.items() if v != None])
    return s

def stringToDict(s):
    """Construct a dictionary from a string"""
    ss = s.split()
    d = dict(list(zip(ss[0::2], ss[1::2])))
    return d


def countLines(file):
    """Count the number of lines in a file"""
    with open(file) as f:
        n = sum(1 for _ in f)
    return n


def directorySetup(config, subdirs=["config", "data", "log", "plots", "report"]):
    for d in subdirs:
        if not os.path.exists(d):
            os.mkdir(d)
        config[d + "_prefix"] = os.path.join(d, config["out_prefix"])

    return config


# cluster configuration is read from json into nested dictionaries
# regular dictionary update loses default values below the first level
# https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
def update(d, u):
    ld = deepcopy(d)
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            if len(v) == 0:
                ld[k] = u[k]
            else:
                r = update(d.get(k, {}), v)
                ld[k] = r
        else:
            ld[k] = u[k]
    return ld


# parent class to represent a compute cluster environment
class Cluster(object):
    """ """
    # constructor
    def __init__(self, std_cluster_file, opt_cluster_file=None, cfg_version="3", verbose=False):
        self.verbose = verbose
        self.class_name = self.__class__.__name__
        # set default pipeline path
        self.pipelinePath = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.submitPath = self.pipelinePath
        # get cluster cfg
        self.openClusterCfg(std_cluster_file, opt_cluster_file, cfg_version, verbose)
        # if resuming, create resume subdirectory
        self.resumeDir = "./resume/"
        if self.clusterCfg["enable_resume"]:
            if not os.path.exists(self.resumeDir):
                os.mkdir(self.resumeDir)

    def analysisInit(self, print_only=False):
        # get analysis name
        self.analysis = os.path.splitext(os.path.basename(os.path.abspath(sys.argv[0])))[0]
        # get user name
        self.username = getpass.getuser()
        # command line
        self.analysisCmd = " ".join(sys.argv[:])
        # get start time (utc)
        tFmt = "%a, %d %b %Y %I:%M:%S %p"
        starttime = datetime.datetime.utcnow()
        self.analysisStart = starttime.strftime(tFmt)
        # for python/shell need "'feb 22 2022'"
        self.analysisStartSec = self.analysisStart.replace(" ", "_")
        # tag for analysis log file name
        self.analysisTag = str(int(time.time()*100))
        # analysis log file
        self.analysisLogFile = "analysis_" + self.analysis + "_" + self.username + \
                                "_" + self.analysisTag + ".log"
        if print_only:
            # print out Info
            print("+++++++++  Print Only +++++++++++")
            print("Analysis: " + self.analysis)
            print("Analysis log file: " + self.analysisLogFile)
            print(self.analysis + " start time: " + self.analysisStart)
        else:
            # open file and output start time
            self.printVerbose("creating analysis log file: " + self.analysisLogFile)
            with open(self.analysisLogFile, "w") as afile:
                afile.write("Analysis: " + self.analysis + "\n")
                afile.write("Cmd: " + self.analysisCmd + "\n")
                afile.write("Version: " + __version__ + "\n")
                afile.write("Start time: " + self.analysisStart + "\n")

    def analysisLog(self, message, print_only = False):
        # append a message to the analysis log file
        if not print_only:
            with open(self.analysisLogFile, "a") as afile:
                afile.write(message + "\n")
        else:
            print(message)


    def getAnalysisName(self):
        return self.analysis

    def getAnalysisLog(self):
        return self.analysisLogFile

    def getAnalysisStart(self):
        return self.analysisStart

    def getAnalysisStartSec(self):
        # get the datetime start (forget the "Sec" part of name)
        return self.analysisStartSec

    def openClusterCfg(self, stdCfgFile, optCfgFile, cfg_version, verbose):
        # get the standard cluster cfg
        self.clusterfile =  os.path.join(self.pipelinePath, stdCfgFile)
        self.printVerbose("reading internal cfg file: " + self.clusterfile)

        with open(self.clusterfile) as cfgFileHandle:
            clusterCfg= json.load(cfgFileHandle)
        # check version
        key = "version"
        if key in clusterCfg:
            if clusterCfg[key] != cfg_version:
                print( "Error: version of : " + stdCfgFile + " should be " + cfg_version +
                       " not " + clusterCfg[key])
                sys.exit(2)
        else:
            print( "Error: version missing in " + stdCfgFile )
            sys.exit(2)
        if self.verbose:
            debugCfg = True
        else:
            key = "debug"
            debugCfg = False
            if key in clusterCfg:
                if clusterCfg[key] == 1:
                    debugCfg = True
        self.clusterCfg = clusterCfg["configuration"]
        if debugCfg:
            print("0>>> Dump of " + clusterCfg["name"] + " ... \n")
            print(json.dumps(self.clusterCfg, indent=3, sort_keys=True))
        if optCfgFile != None:
            self.printVerbose("reading user cfg file: " + optCfgFile)

            with open(optCfgFile) as cfgFileHandle:
                clusterCfg = json.load(cfgFileHandle)
            optCfg = clusterCfg["configuration"]
            if debugCfg:
                print("0>>> Dump of opt cfg file " + optCfgFile + " ... \n")
                print(json.dumps(optCfg, indent=3, sort_keys=True))
            # update
            self.clusterCfg = update(self.clusterCfg, optCfg)
            if debugCfg:
                print("0>>> Dump of updated cluster cfg ... \n")
                print(json.dumps(self.clusterCfg, indent=3, sort_keys=True))
        key = "memory_limits"
        if key not in self.clusterCfg:
            self.clusterCfg[key] = None

        # update pipeline path if specified
        key = "pipeline_path"
        if key in self.clusterCfg:
            self.pipelinePath = self.clusterCfg["pipeline_path"]

    def getPipelinePath(self):
        return self.pipelinePath

    def getSubmitPath(self):
        return self.submitPath

    def getClusterCfg(self):
        return self.clusterCfg

    def memoryLimit(self, job_name):
        memlim = None
        memLimits = self.clusterCfg["memory_limits"]
        if memLimits is None:
            return memlim
        # find a dicitionary of all partial match with job_name
        fd = dict([elem for elem in list(memLimits.items()) if job_name.find(elem[0]) !=-1])
        # if we have an exact match
        if job_name in list(fd.keys()):
            memlim = fd[job_name]
        # else we have a partial match
        else:
            jobMem = [ v for k,v in fd.items() if job_name.find(k) != -1 ]
            if len(jobMem):
                # find if there is an exact match with jobname and memlimit's key
                # just find the first match to job_name
                memlim = jobMem[0]
        self.printVerbose('Memory Limit - job: ' + job_name + " memlim: " +
                          str(memlim) + "MB")
        return memlim

    def printVerbose(self, message):
        if self.verbose:
            print(">>> " + self.class_name + ": " + message)


class Docker_Cluster(Cluster):
    def __init__(self, std_cluster_file, opt_cluster_file=None, cfg_version="3", verbose=True):
        self.class_name = self.__class__.__name__
        self.std_cluster_file = std_cluster_file
        super(Docker_Cluster, self).__init__(std_cluster_file, opt_cluster_file, cfg_version, verbose)

    def runCmd(self, job_name, cmd, logfile=None):
        # get the docker image
        rc = self.clusterCfg["run_cmd"]
        dockerimage = rc["docker_image"]
        # build the docker run command to run the R driver and all it's args
        working_dir = os.getenv('PWD')
        drun = "docker run --rm -it -v /projects:/projects -w " + working_dir + " " + \
                dockerimage + " "
        dargs = cmd[0] + ' ' + " ".join(cmd[1:])
        dcmd = drun + dargs

        self.printVerbose("docker run cmd:\n\t" + dcmd)
        (pmsg, status) = port_popen.popen_stdout(dcmd, logfile, jobname="run_cmd")
        if status != 0:
            print(">>> Error: Docker_Cluster/runCmd executing popen\nStatus: " + str(status) + " Msg: " + pmsg)
            sys.exit(status)


class AWS_Batch(Docker_Cluster):

    def __init__(self, opt_cluster_file=None, verbose=False):
        if not batchSupport:
            print('AWS_Batch: AWS Batch is not supported (boto3 is not available)')
            sys.exit(2)
        self.class_name = self.__class__.__name__
        self.std_cluster_file = "./aws_batch_cfg.json"
        self.opt_cluster_file = opt_cluster_file
        cfgVersion = "3.3"
        super(AWS_Batch, self).__init__(self.std_cluster_file, opt_cluster_file, cfgVersion, verbose)

        # get the job parameters
        self.jobParams = self.clusterCfg["job_parameters"]
        #user = getpass.getuser()
        wdkey = "wd"
        if wdkey not in self.jobParams or self.jobParams[wdkey] == "":
            self.jobParams[wdkey] = os.getenv('PWD')

        # set maxperf
        self.maxperf = True
        if self.clusterCfg["maxperf"] == 0:
            self.maxperf = False

        # get the submit options
        self.submitOpts = self.clusterCfg["submit_opts"]

        # get the run cmd options
        self.runCmdOpts = self.clusterCfg["run_cmd"]

        # get the sync job options
        self.syncOpts = self.clusterCfg["sync_job"]

        # get the queue
        self.queue = self.clusterCfg["queue"]

        # create the batch client
        try:
            session = boto3.Session(profile_name = self.clusterCfg["aws_profile"])
            self.batchC = session.client('batch')
        except Exception as e:
            print('boto3 session or client exception ' + str(e))
            sys.exit(2)

        # retryStrategy
        self.retryStrategy = self.clusterCfg["retryStrategy"]

    def analysisInit(self, print_only=False):
        # base init first
        super(AWS_Batch, self).analysisInit(print_only)
        # jobinfo file name
        self.jiFileName = self.analysis + "_" + self.analysisTag + "_jobinfo.txt"
        # analysis log file
        profile = self.clusterCfg["aws_profile"]
        if print_only:
            print("+++++++++  Print Only +++++++++++")
            print("AWS profile: " + profile)
            print("AWS job def: " + self.submitOpts["jobdef"])
            print("AWS batch queue: " + self.queue)
        else:
            self.printVerbose("AWS profile: " + profile)
            self.printVerbose("AWS job def: " + self.submitOpts["jobdef"])
            self.printVerbose("AWS batch queue: " + self.queue)
            with open(self.analysisLogFile, "a") as afile:
                afile.write("AWS profile: " + profile + "\n")
                afile.write("AWS job def: " + self.submitOpts["jobdef"] + "\n")
                afile.write("AWS batch queue: " + self.queue + "\n")

    def submitJob(self, job_name, cmd, args=None, holdid=None, array_range=None,
                  request_cores=None, print_only=False, **kwargs):
        # change args from list to str
        if args is None:
            args = ["NoArgs"]
        args_str = " ".join(args)
        # fill the submit parameters dict
        awsbatch.subParams['profile'] = self.clusterCfg["aws_profile"]
        awsbatch.subParams['cluster_file']= self.opt_cluster_file
        awsbatch.subParams['clustercfg'] = self.clusterCfg
        awsbatch.subParams['jobname'] = job_name
        awsbatch.subParams['cmd'] = os.path.basename(cmd)
        awsbatch.subParams['args'] = args_str
        awsbatch.subParams['holdid'] = holdid
        awsbatch.subParams['array_range'] = array_range
        awsbatch.subParams['request_cores'] = request_cores
        awsbatch.subParams['print_only'] = print_only
        awsbatch.subParams['verbose'] = self.verbose
        awsbatch.subParams['maxmem'] = None
        awsbatch.subParams['workdir'] = None
        awsbatch.subParams['queue'] = None
        awsbatch.subParams['jobdef'] = None
        awsbatch.subParams['apath'] = os.path.dirname(cmd)
        awsbatch.subParams['infofile'] = self.jiFileName
        awsbatch.subParams['analysislog'] = "update"
        # submit the job
        submit_id = awsbatch.submitjob(awsbatch.subParams)
        # update analysis log
        super(AWS_Batch, self).analysisLog(awsbatch.subParams['analysislog'], print_only)

        return submit_id

class SGE_Cluster(Cluster):

    def __init__(self, std_cluster_file, opt_cluster_file=None, cfg_version="3", verbose=True):
        self.class_name = self.__class__.__name__
        self.std_cluster_file = std_cluster_file
        super(SGE_Cluster, self).__init__(std_cluster_file, opt_cluster_file, cfg_version, verbose)

    def runCmd(self, job_name, cmd, logfile=None):
        # get and set the env
        key = "-v"
        if key in list(self.clusterCfg["submit_opts"].keys()):
            vopt = self.clusterCfg["submit_opts"][key]
            envVars = vopt.split(",")
            for var in envVars:
                varVal = var.split("=")
                check = "$PATH"
                if varVal[1].endswith(check):
                    np = varVal[1][:-len(check)]
                    cp = os.environ['PATH']
                    os.environ[varVal[0]] = np + cp
                else:
                    os.environ[varVal[0]] = varVal[1]

        (pmsg, status) = port_popen.popen_stdout(cmd, logfile, jobname=job_name, shell=False)
        if status != 0:
            print(">>> Error: SGE_Cluster:runCmd executing popen\nStatus: " + str(status) + " Msg: " + str(pmsg))
            sys.exit(status)

    def submitJob(self, binary=False, hold_array=None, **kwargs):
        subOpts = deepcopy(self.clusterCfg["submit_opts"])
        # set the submit cmd (e.g., qsub)
        submit_cmd = self.clusterCfg["submit_cmd"]
        # check for binary
        if binary:
            subOpts["-b"] = "y"
        # job name
        job_name = kwargs["job_name"]
        subOpts["-N"] = job_name
        lmsg = "Job: " + job_name
        # if not holding for job array (element wise), check for holding for jobs
        if hold_array == None:
            key = "holdid"
            if key in kwargs and kwargs[key] != []:
                if isinstance(kwargs[key], str):
                    kwargs[key] = [kwargs[key]]
                subOpts["-hold_jid"] =  ",".join(kwargs[key])
        else:
            subOpts["-hold_jid_ad"] = hold_array
        # array job
        key = "array_range"
        array_job = False
        lmsg_array = "no"
        if key in kwargs:
            subOpts["-t"] = kwargs[key]
            lmsg_array = kwargs[key]
            array_job = True
        lmsg = lmsg + " /array: " + lmsg_array
        # threads
        key = "request_cores"
        lmsg_cores = "1"
        memcoreFactor = 1.
        if key in kwargs and kwargs[key] != None and kwargs[key] != "1":
            # adjust memcoreFactor if a specific no. of cores is passed (i.e., no "-")
            reqCores = kwargs[key]
            if not "-" in reqCores:
                memcoreFactor = float(reqCores)
            subOpts["-pe"] = self.clusterCfg["parallel_env"] + " " + reqCores
            lmsg_cores = reqCores
        lmsg = lmsg + " /cores: " + lmsg_cores
        # get memory limit option (adjust based on specifying a specific number of cores)
        key = "memory_limits"
        lmsg_mem = "not provided"
        if key in list(self.clusterCfg.keys()):
            memlim = super(SGE_Cluster, self).memoryLimit(job_name)
            if memlim != None:
                memlim = memlim/memcoreFactor
                subOpts["-l"] = "h_vmem="+str(memlim)+"M"
                lmsg_mem = str(memlim)
        lmsg = lmsg + " /memlim: " + lmsg_mem
        # email
        key = "email"
        if key in kwargs and kwargs[key] != None:
            subOpts["-m"] = "e"
            subOpts["-M"] = kwargs[key]
        # driver (cmd) args
        key = "args"
        if not key in kwargs:
            kwargs[key] = []
        argStr = " ".join(kwargs[key])
        # check for eqw
        if self.clusterCfg["enable_eqw"]:
            eqwVal = "ENABLE_EQW=1"
            key = "-v"
            if key in subOpts:
                subOpts[key] = subOpts[key] + "," + eqwVal
            else:
                subOpts[key] = eqwVal
        # create a string for the submit options
        optStr = dictToString(subOpts)
        # create the entire submit command
        sub_cmd = " ".join([submit_cmd, optStr, kwargs["cmd"], argStr])
        # check for resume
        if self.clusterCfg["enable_resume"]:
            rscript = self.clusterCfg["resume_script"]
            spath = super(SGE_Cluster, self).getSubmitPath()
            rscript = spath + "/" + rscript + " sge " + job_name
            # update submit opts -N
            subOpts["-N"] = "resume_" + subOpts["-N"]
            subOpts["-b"] = "y"
            subOpts["-o"] = self.resumeDir
            optStr = dictToString(subOpts)
            # update sub_cmd
            sub_cmd = " ".join([submit_cmd, optStr, rscript, kwargs["cmd"], argStr])

        key = "print_only"
        if key in kwargs and kwargs[key] == True:
            print(sub_cmd)
            return "000000"
        self.printVerbose("submitting job: " + sub_cmd)
        (jobid, status) = port_popen.popen(sub_cmd)
        if status != 0:
            print(">>> Error: SGE_Cluster:submitJob executing popen\nStatus: " + str(status) + " Msg: " + str(jobid))
            sys.exit(status)

        if array_job:
            jobid = jobid.split(".")[0]
        print("Submitting job " + jobid + " (" + subOpts["-N"] + ")")
        lmsg += " /jobid: " + str(jobid)
        super(SGE_Cluster, self).analysisLog(lmsg)

        return jobid

class UW_Cluster(SGE_Cluster):

    def __init__(self, opt_cluster_file=None, verbose=False):
        self.class_name = self.__class__.__name__
        self.std_cluster_file = "./cluster_cfg.json"
        cfgVersion="3"
        super(UW_Cluster, self).__init__(self.std_cluster_file, opt_cluster_file, cfgVersion, verbose)


class AWS_Cluster(SGE_Cluster):

    def __init__(self, opt_cluster_file=None, verbose=False):
        self.class_name = self.__class__.__name__
        self.std_cluster_file = "./aws_cluster_cfg.json"
        cfgVersion="3"
        super(AWS_Cluster, self).__init__(self.std_cluster_file, opt_cluster_file, cfgVersion, verbose)

    def submitJob(self, **kwargs):
        # currently, no email on aws
        kwargs["email"] = None
        jobid = super(AWS_Cluster, self).submitJob(**kwargs)
        return jobid

class Slurm_Cluster(Docker_Cluster):

    def __init__(self, opt_cluster_file=None, verbose=False):
        self.class_name = self.__class__.__name__
        self.std_cluster_file = "./slurm_cfg.json"
        self.opt_cluster_file = opt_cluster_file

        cfgVersion = "2.2"
        super(Slurm_Cluster, self).__init__(self.std_cluster_file, opt_cluster_file, cfgVersion, verbose)
        # open slurm partitions cfg
        self.openPartitionCfg(self.pipelinePath + "/" + self.clusterCfg["partition_cfg"])
        # update pipelinePath
        key = "pipeline_path_docker"
        if key in list(self.clusterCfg.keys()):
            self.pipelinePath = self.clusterCfg[key]

    def openPartitionCfg(self, a_pcfg):
        # open partition cfg and set partition names, partitions for cluster
        with open(a_pcfg) as cfgFileHandle:
            partitionCfg = json.load(cfgFileHandle)
        cname = self.clusterCfg["cluster"]
        if cname not in list(partitionCfg["clusters"].keys()):
            print("Cluster " + cname + " not found in " + a_pcfg)
            sys.exit(2)
        self.partitions = partitionCfg["clusters"][cname]
        self.partition_names = list(self.partitions.keys())

    def analysisInit(self, print_only=False):
        # analysis log file and analysis info
        super(Slurm_Cluster, self).analysisInit(print_only)
        super(Slurm_Cluster, self).analysisLog("Slurm cluster: " +
                                               self.clusterCfg["cluster"], print_only)
        super(Slurm_Cluster, self).analysisLog("Slurm submit script: " +
                                               self.clusterCfg["submit_script"] + "\n", print_only)

    def getPartition(self, a_jobname, a_memsize, a_reqcores, a_tasksPerPartition):
        # find all partitions with mem >= a_memsize*tpp
        memcheck = a_memsize*a_tasksPerPartition
        pmem = [ k for k in self.partition_names if self.partitions[k]["mem"] > memcheck ]
        if len(pmem) == 0:
            print("Error: cannot find partition with sufficient memory (" + str(memcheck) + "MB)")
            sys.exit(2)
        # from partitions with mem, find partitions with cores >= a_reqcores*tpp
        corecheck = a_reqcores*a_tasksPerPartition
        pmemcore = [ k for k in pmem if self.partitions[k]["cores"] >= corecheck ]
        if len(pmemcore) == 0:
            print("Error: cannot find partition with enough cores (" + str(corecheck) +
                 ") for "+ str(memcheck) + "MB memory")
            sys.exit(2)

        # from partitions with mem & core, find partition with min memory
        thepart = pmemcore[0]
        themem = self.partitions[thepart]["mem"]
        nop = len(pmemcore)
        if nop > 1 :
            for i in range(1,nop):
                pcheck = pmemcore[i]
                mcheck = self.partitions[pcheck]["mem"]
                if mcheck < themem:
                    thepart = pcheck
                    themem = mcheck
        return thepart

    def submitJob(self, **kwargs):
        # get the various config attributes
        submitOpts = deepcopy(self.clusterCfg["submit_opts"])
        dockerOpts = deepcopy(self.clusterCfg["rdocker_opts"])
        cluster = self.clusterCfg["cluster"]
        job_cmd = self.clusterCfg["submit_cmd"]
        dependency_value = self.clusterCfg["dependency_value"]
        pipeline_path_docker = self.clusterCfg["pipeline_path_docker"]
        # get tasks per partition based on job name
        tasksPerPartition = 1;
        tppDict = self.clusterCfg["tasks_per_partition"]
        jobName = kwargs["job_name"]
        jobPart = [ v for k,v in tppDict.items() if jobName.find(k) != -1]
        if len(jobPart):
            # just find the first match to jobname
            tasksPerPartition = jobPart[0]
        self.printVerbose("tasks per partition for job " + jobName + ": " +
                          str(tasksPerPartition))
        # process kwargs for submit options
        submitOpts["--job-name"] = jobName
        lmsg = "Job: " + jobName

        key = "holdid"
        if key in kwargs and kwargs[key] != []:
            if isinstance(kwargs[key], str):
                kwargs[key] = [kwargs[key]]
            submitOpts["--dependency"] =  dependency_value + ":" + ":".join(kwargs[key])

        key = "array_range"
        lmsg_array = "no"
        if key in kwargs:
            submitOpts["--array"] = kwargs[key]
            lmsg_array = kwargs[key]
        lmsg = lmsg + " /array: " + lmsg_array

        key = "request_cores"
        lmsg_cores = "1"
        reqCores = 1
        if key in kwargs and kwargs[key] != None and kwargs[key] != "1":
            # just take the max if specified 1-max
            rcl = kwargs[key].split("-")
            reqCores = rcl[len(rcl)-1]
            submitOpts["--cpus-per-task"] = reqCores
            lmsg_cores = reqCores
        lmsg = lmsg + " /cores: " + lmsg_cores

        # get memory limit option
        key = "memory_limits"
        lmsg_mem = "not provided"
        if key in list(self.clusterCfg.keys()):
            memlim = super(Slurm_Cluster, self).memoryLimit(kwargs["job_name"])
            if memlim == None:
                memlim = 8000
            submitOpts["--mem"] = str(memlim) + "M"
            lmsg_mem = submitOpts["--mem"]
        lmsg = lmsg + " /memlim: " + lmsg_mem
        dockerOpts["--mem_limit"] = str(int(math.ceil(float(memlim)/1000)))

        # get partition
        thePartition = self.getPartition(kwargs["job_name"], memlim, int(reqCores), tasksPerPartition)
        submitOpts["--partition"] = thePartition
        lmsg = lmsg + " /cluster: " + cluster
        lmsg = lmsg + " /parition: " + submitOpts["--partition"]
        lmsg = lmsg + " /tasks_per_partition: " + str(tasksPerPartition)
        # get the machine and cost
        theMachine = self.partitions[thePartition]["machine"]
        theCost = str((self.partitions[thePartition]["cost"]/tasksPerPartition))
        lmsg = lmsg + "/machine: " + theMachine + " ( " + theCost + "/hr for " + str(tasksPerPartition) + " task(s))"
        # submit output (log)
        submit_logdir = ""
        if self.clusterCfg["submit_log_dir"] != None:
            submit_logdir = self.clusterCfg["submit_log_dir"] + "/sbatch_"
        if submitOpts["--array"] == None:
            submitOpts["--output"] = submit_logdir + submitOpts["--job-name"] + "_%j.log"
        else:
            submitOpts["--output"] = submit_logdir + submitOpts["--job-name"] + "_%A_%a.log"

        # set the docker options
        # -- cmd (change path associated within docker)
        cmd = kwargs["cmd"]
        dockerOpts["--runcmd"] = pipeline_path_docker + "/" + os.path.basename(cmd)
        # -- args
        key = "args"
        if not key in kwargs:
            kwargs[key] = []
        dockerOpts["--runargs"] = '"' + " ".join(kwargs[key]) + '"'
        # -- cost
        dockerOpts["--machine"] = theMachine
        dockerOpts["--cost"] = theCost
        # -- working_dir
        dockerOpts["--working_dir"] = os.getenv('PWD')
        # dockerOpts - boolean options.  All boolean options are default to false
        # in runDocker.py.  A boolean option does not have a value.
        # If a boolean option is passed then that option in runDocker is set true.
        # We convert all options to strings so to a true boolean is passed as an empty
        # string; for false option, we delete the key.
        # -- verbose
        key = "--verbose"
        if self.verbose:
            dockerOpts[key] = True
        if dockerOpts[key]:
            dockerOpts[key] = ""
        else:
            del dockerOpts[key]
        # -- stats
        key = "--stats"
        if dockerOpts[key]:
            dockerOpts[key] = ""
        else:
            del dockerOpts[key]
        # -- pull
        key = "--pull"
        if dockerOpts[key]:
            dockerOpts[key] = ""
        else:
            del dockerOpts[key]
        # -- log
        key = "--log"
        if dockerOpts[key]:
            dockerOpts[key] = ""
        else:
            del dockerOpts[key]
        # convert opts to strings
        suboptStr = dictToString(submitOpts)
        dockeroptStr = dictToString(dockerOpts)
        # submit scripts
        submit_script = self.submitPath + "/" + self.clusterCfg["submit_script"]
        submit_prescript = self.submitPath + "/" + self.clusterCfg["submit_prescript"]
        resume_script = self.submitPath + "/" + self.clusterCfg["resume_script"]
        if self.clusterCfg["enable_resume"]:
            # resume
            submit_script = submit_prescript + " " + resume_script + " slurm " + jobName + " " + submit_script
            # update submit opts - batch job name and log file
            rpre = "resume_"
            submitOpts["--job-name"] = rpre + jobName
            # either a single job or an array job
            if submitOpts["--array"] == None:
                submitOpts["--output"] = self.resumeDir + submitOpts["--job-name"] + "_%j.log"
            else:
                submitOpts["--output"] = self.resumeDir + submitOpts["--job-name"] + "_%A_%a.log"

            suboptStr = dictToString(submitOpts)
        else:
            self.clusterCfg["execute_prescript"]
            submit_script = submit_prescript + " " + submit_script

        sub_cmd = " ".join([job_cmd, suboptStr, submit_script, dockeroptStr])

        key = "print_only"
        po = False
        if key in kwargs and kwargs[key] == True:
            po = True

        super(Slurm_Cluster, self).analysisLog(lmsg, po)
        if po:
            print(sub_cmd + "\n")
            jobid = submitOpts["--job-name"]
        else:
            self.printVerbose("submitting job: " + sub_cmd)
            super(Slurm_Cluster, self).analysisLog("> sbatch: " + sub_cmd)
            (sub_out, status) = port_popen.popen(sub_cmd)
            if status != 0:
                print(">>> Error: Slurm_Cluster:submitJob executing popen\nStatus: " + \
                      str(status) + " Msg: " + str(sub_out))
                sys.exit(status)
            jobid = sub_out.split(" ")[3].strip()
            super(Slurm_Cluster, self).analysisLog("> jobid: " + str(jobid) + "\n")
            print(sub_out + " cluster: " + cluster + " / job: " + submitOpts["--job-name"] +
                  " / job id: " + jobid)
        return jobid

class GCP_Cluster(Slurm_Cluster):

    def __init__(self, opt_cluster_file=None, verbose=False):
        self.class_name = self.__class__.__name__
        self.std_cluster_file = "./gcp_slurm_cfg.json"
        cfgVersion="3"
        super(GCP_Cluster, self).__init__(self.std_cluster_file, opt_cluster_file, cfgVersion, verbose)

class ClusterFactory(object):

    @staticmethod
    def createCluster(cluster_type, cluster_file=None, verbose=False):
        allSubClasses = getAllSubclasses(Cluster)
        for subclass in allSubClasses:
            if subclass.__name__ == cluster_type:
                return subclass(cluster_file, verbose)
        raise Exception("unknown cluster type: " + cluster_type + "!")

def getAllSubclasses(base):
    all_subclasses = []
    for subclass in base.__subclasses__():
        all_subclasses.append(subclass)
        all_subclasses.extend(getAllSubclasses(subclass))
    return all_subclasses
