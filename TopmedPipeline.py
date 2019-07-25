"""Utility functions for TOPMed pipeline"""

__version__ = "2.4.3"

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
from datetime import datetime, timedelta
import awsbatch

try:
    import boto3
    import batchInit
except ImportError:
    print ("AWS batch not supported.")



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
    for key, value in config.iteritems():
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
        dummy = reader.next()
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
    return range(start, end + 1)

def parseChromosomes(chromosomes):
    chromString = " ".join([str(x) for x in chromosomeRangeToList(chromosomes)])
    chromString = chromString.replace("23", "X")
    chromString = chromString.replace("24", "Y")
    return chromString


def dictToString(d):
    """Construct a string from a dictionary"""
    s = ' '.join([k + ' ' + v for k, v in d.iteritems()])
    return s

def stringToDict(s):
    """Construct a dictionary from a string"""
    ss = s.split()
    d = dict(zip(ss[0::2], ss[1::2]))
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
    for k, v in u.iteritems():
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
        # set analysis name

        self.openClusterCfg(std_cluster_file, opt_cluster_file, cfg_version, verbose)

    def analysisInit(self, print_only=False):
        # get analysis name
        self.analysis = os.path.splitext(os.path.basename(os.path.abspath(sys.argv[0])))[0]
        # get user name
        self.username = getpass.getuser()
        # ascii time
        self.analysisStart = time.asctime()
        # command line
        self.analysisCmd = " ".join(sys.argv[:])
        # sec time
        dt_ref = datetime(1970,1,1)
        tFmt = '%a %b %d %H:%M:%S %Y'
        dt_start = datetime.strptime(self.analysisStart, tFmt)
        self.analysisStartSec = str((dt_start - dt_ref).total_seconds())
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
            self.printVerbose("Creating analysis log file: " + self.analysisLogFile)
            with open(self.analysisLogFile, "w") as afile:
                afile.write("Analysis: " + self.analysis + "\n")
                afile.write("Cmd: " + self.analysisCmd + "\n")
                afile.write("Version: " + __version__ + "\n")
                afile.write("Start time: " + self.analysisStart + "\n")

    def analysisLog(self, message):
        # append a message to the analysis log file
        with open(self.analysisLogFile, "a") as afile:
            afile.write("Message: " + message + "\n")

    def getAnalysisName(self):
        return self.analysis

    def getAnalysisLog(self):
        return self.analysisLogFile

    def getAnalysisStart(self):
        return self.analysisStart

    def getAnalysisStartSec(self):
        return self.analysisStartSec

    def openClusterCfg(self, stdCfgFile, optCfgFile, cfg_version, verbose):
        # get the standard cluster cfg
        self.clusterfile =  os.path.join(self.pipelinePath, stdCfgFile)
        self.printVerbose("0>> openClusterCfg: Reading internal cfg file: " + self.clusterfile)

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
            print json.dumps(self.clusterCfg, indent=3, sort_keys=True)
        if optCfgFile != None:
            self.printVerbose("0>> openClusterCfg: Reading user cfg file: " + optCfgFile)

            with open(optCfgFile) as cfgFileHandle:
                clusterCfg = json.load(cfgFileHandle)
            optCfg = clusterCfg["configuration"]
            if debugCfg:
                print("0>>> Dump of " + clusterCfg["name"] + " ... \n")
                print json.dumps(optCfg, indent=3, sort_keys=True)
            # update
            self.clusterCfg = update(self.clusterCfg, optCfg)
            if debugCfg:
                print("0>>> Dump of updated cluster cfg ... \n")
                print json.dumps(self.clusterCfg, indent=3, sort_keys=True)
        key = "memory_limits"
        if key not in self.clusterCfg:
            self.clusterCfg[key] = None

        # update pipeline path if specified
        key = "pipeline_path"
        if key in self.clusterCfg:
            self.pipelinePath = self.clusterCfg["pipeline_path"]

    def getPipelinePath(self):
        return self.pipelinePath

    def getClusterCfg(self):
        return self.clusterCfg

    def memoryLimit(self, job_name):
        memlim = None
        memLimits = self.clusterCfg["memory_limits"]
        if memLimits is None:
            return memlim
        # find a dicitionary of all partial match with job_name
        fd = dict(filter(lambda elem: job_name.find(elem[0]) !=-1,memLimits.items()))
        # if we have an exact match
        if job_name in fd.keys():
            memlim = fd[job_name]
        # else we have a partial match
        else:
            jobMem = [ v for k,v in fd.iteritems() if job_name.find(k) != -1 ]
            if len(jobMem):
                # find if there is an exact match with jobname and memlimit's key
                # just find the first match to job_name
                memlim = jobMem[0]
        self.printVerbose('\t>>> Memory Limit - job: ' + job_name + " memlim: " +
                          str(memlim) + "MB")
        return memlim

    def printVerbose(self, message):
        if self.verbose:
            print(message)


class AWS_Batch(Cluster):

    def __init__(self, opt_cluster_file=None, verbose=False):
        self.class_name = self.__class__.__name__
        self.std_cluster_file = "./aws_batch_cfg.json"
        self.opt_cluster_file = opt_cluster_file
        cfgVersion = "3.2"
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
        # analysis log file and autogen stuff
        profile = self.clusterCfg["aws_profile"]
        # autogen (for auto generation of queue and ce)
        self.autogen_ce = False
        if self.clusterCfg["autogen_ce"].lower() == "yes":
            self.autogen_ce = True
        if self.autogen_ce:
            # pricing
            pricing = self.clusterCfg["pricing"]
            # create an aws tag for costs
            self.awstag = self.analysis + "_" + self.username + "_" + self.analysisTag
            # create ce name
            ceName = "ag_ce" + "_" + pricing + "_" + self.analysis + "_" + self.username
            # create costing queue name
            awsqueue = "ag_queue" + "_" + pricing + "_" + self.analysis + "_" + self.username
            self.queue = awsqueue
            print("Creating batch env ...")
            # if print_only
            if print_only:
                print("+++++++++  Print Only +++++++++++")
                print("AWS autogen: " + str(self.autogen_ce))
                print("AWS tag: " + self.awstag)
                print("AWS profile: " + profile)
                print("AWS batch queue: " + self.queue)
                print("AWS batch ce: " + ceName)
                print("AWS job def: " + self.submitOpts["jobdef"])
                print("Check AWS batch queue to verify ...")
            else:
                with open(self.analysisLogFile, "a") as afile:
                    afile.write("AWS profile: " + profile + "\n")
                    afile.write("AWS batch queue: " + self.queue + "\n")
                    afile.write("AWS batch ce: " + ceName + "\n")
                    afile.write("AWS tag: " + self.awstag + "\n")
            # set default instance types
            instanceTypes = None
            if len(self.clusterCfg["batch_instances"]) > 0:
                instanceTypes = self.clusterCfg["batch_instances"]
            # init batch by creating ce and queue based on batch pricing
            if not print_only:
                self.printVerbose("AWS tag: " + self.awstag)
                self.printVerbose("AWS profile: " + profile)
                self.printVerbose("batch queue: " + self.queue)
                self.printVerbose("batch pricing: " + self.clusterCfg["pricing"])
                self.printVerbose("AWS job def: " + self.submitOpts["jobdef"])
                self.printVerbose("instance types: " + str(instanceTypes))
                self.printVerbose("compute environment name: " + ceName)
                self.printVerbose("aws profile: " + profile)
            batchInit.createEnv(self.batchC,
                                self.queue,
                                self.clusterCfg["pricing"],
                                instanceTypes,
                                ceName,
                                self.awstag,
                                profile,
                                self.verbose)
        else:
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

    def runCmd(self, job_name, cmd, logfile=None):
        # redirect stdout/stderr
        self.printVerbose("1===== runCmd: job " + job_name + " beginning ...")
        self.printVerbose("1===== runCmd: cmd " + str(cmd))
        if logfile != None:
            sout = sys.stdout
            serr = sys.stderr
            flog = open ( logfile, 'w' )
            sys.stderr = sys.stdout = flog
        # spawn
        self.printVerbose("1===== runCmd: executing " + str(cmd))
        process = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr, shell=False)
        status = process.wait()
        # redirect stdout back
        if logfile != "":
            sys.stdout = sout
            sys.stderr = serr
        if status:
            print( "Error running job: " + job_name + " (" + str(status) + ") for command:" )
            print( "\t> " + str(cmd) )
            sys.exit(2)

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
        super(AWS_Batch, self).analysisLog(awsbatch.subParams['analysislog'])

        return submit_id

class SGE_Cluster(Cluster):

    def __init__(self, std_cluster_file, opt_cluster_file=None, cfg_version="3", verbose=True):
        self.class_name = self.__class__.__name__
        self.std_cluster_file = std_cluster_file
        super(SGE_Cluster, self).__init__(std_cluster_file, opt_cluster_file, cfg_version, verbose)

    def runCmd(self, job_name, cmd, logfile=None):
        # get and set the env
        key = "-v"
        if key in self.clusterCfg["submit_opts"].keys():
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
        # redirect stdout/stderr
        if logfile != None:
            sout = sys.stdout
            serr = sys.stderr
            flog = open ( logfile, 'w' )
            sys.stderr = sys.stdout = flog
        # spawn
        process = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr, shell=False)
        status = process.wait()
        # redirect stdout back
        if logfile != "":
            sys.stdout = sout
            sys.stderr = serr
        if status:
            print( "Error running job: " + job_name + " (" + str(status) + ") for command:" )
            print( "\t> " + str(cmd) )
            sys.exit(2)

    def submitJob(self, **kwargs):
        subOpts = deepcopy(self.clusterCfg["submit_opts"])
        # set the job cmd
        kwargs["job_cmd"] = self.clusterCfg["submit_cmd"]

        jobid = self.executeJobCmd(subOpts, **kwargs)
        return jobid

    def executeJobCmd(self, submitOpts, **kwargs):
        # update sge opts
        submitOpts["-N"] = kwargs["job_name"]
        lmsg = "Job: " + kwargs["job_name"]

        key = "holdid"
        if key in kwargs and kwargs["holdid"] != []:
            if isinstance(kwargs["holdid"], str):
                kwargs["holdid"] = [kwargs["holdid"]]
            submitOpts["-hold_jid"] =  ",".join(kwargs["holdid"])

        key = "array_range"
        lmsg_array = "no"
        if key in kwargs:
            submitOpts["-t"] = kwargs["array_range"]
            lmsg_array = kwargs["array_range"]
        lmsg = lmsg + " / array: " + lmsg_array

        key = "request_cores"
        lmsg_cores = "1"
        memcoreFactor = 1.
        if key in kwargs and kwargs[key] != None and kwargs[key] != "1":
            # adjust memcoreFactor if a specific no. of cores is passed (i.e., no "-")
            reqCores = kwargs["request_cores"]
            if not "-" in reqCores:
                memcoreFactor = float(reqCores)
            submitOpts["-pe"] = self.clusterCfg["parallel_env"] + " " + reqCores
            lmsg_cores = reqCores
        lmsg = lmsg + " / cores: " + lmsg_cores

        # get memory limit option (adjust based on specifying a specific number of cores)
        key = "memory_limits"
        lmsg_mem = "not provided"
        if key in self.clusterCfg.keys():
            memlim = super(SGE_Cluster, self).memoryLimit(kwargs["job_name"])
            if memlim != None:
                memlim = memlim/memcoreFactor
                submitOpts["-l"] = "h_vmem="+str(memlim)+"M"
                lmsg_mem = str(memlim)
        lmsg = lmsg + " / memlim: " + lmsg_mem


        key = "email"
        if key in kwargs and kwargs[key] != None:
            submitOpts["-m"] = "e"
            submitOpts["-M"] = kwargs["email"]

        # update sge cmd and args
        key = "args"
        if not key in kwargs:
            kwargs["args"] = []
        argStr = " ".join(kwargs["args"])

        optStr = dictToString(submitOpts)

        sub_cmd = " ".join([kwargs["job_cmd"], optStr, kwargs["cmd"], argStr])

        key = "print_only"
        if key in kwargs and kwargs[key] == True:
            print sub_cmd
            return "000000"
        self.printVerbose("executeJobCmd subprocess/sub_cmd " + sub_cmd)
        super(SGE_Cluster, self).analysisLog(lmsg)
        process = subprocess.Popen(sub_cmd, shell=True, stdout=subprocess.PIPE)
        pipe = process.stdout
        sub_out = pipe.readline()
        jobid = sub_out.strip(' \t\n\r')

        if "array_range" in kwargs:
            jobid = jobid.split(".")[0]
        print("Submitting job " + jobid + " (" + kwargs["job_name"] + ")")

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
