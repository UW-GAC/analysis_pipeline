"""Utility functions for TOPMed pipeline"""

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

try:
    import boto3
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
        if line[0][0] == "#":
            continue

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
    for k, v in u.iteritems():
        if isinstance(v, collections.Mapping):
            r = update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d


# parent class to represent a compute cluster environment
class Cluster(object):
    """ """
    # constructor
    def __init__(self, cluster_file=None, verbose=False):
        self.verbose = verbose
        self.class_name = self.__class__.__name__
        # get the standard cluster cfg
        pipePath = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.clusterfile =  os.path.join(pipePath, "./cluster_cfg.json")
        self.printVerbose("0>> ClusterCtor: Reading internal cfg file: " + self.clusterfile)

        with open(self.clusterfile) as cfgFileHandle:
            clusterCfg= json.load(cfgFileHandle)
        self.clusterCfg = clusterCfg["cluster_types"][self.class_name]
        if cluster_file != None:
            self.printVerbose("0>> ClusterCtor: Reading user cfg file: " + self.clusterfile)

            with open(cluster_file) as cfgFileHandle:
                clusterCfg = json.load(cfgFileHandle)
            optCfg = clusterCfg["cluster_types"][self.class_name]
            # update
            self.clusterCfg = update(self.clusterCfg, optCfg)

    def getClusterCfg(self):
        return self.clusterCfg

    def memoryLimit(self, memLimits, job_name):
        memlim = None
        a =[ akey for akey in memLimits.keys() if job_name.startswith( akey ) ]
        if len(a):
            # just find the first match to job_name
            memlim = memLimits[a[0]]
        return memlim

    def printVerbose(self, message):
        if self.verbose:
            print(message)


class AWS_Batch(Cluster):

    def __init__(self, cluster_file=None, verbose=False):
        self.class_name = self.__class__.__name__
        super(AWS_Batch, self).__init__(cluster_file, verbose)

        # get the job parameters
        self.jobParams = self.clusterCfg["job_parameters"]
        # set the working directory to cwd if not in clusterCfg

        #user = getpass.getuser()
        wdkey = "wd"
        if wdkey not in self.jobParams or self.jobParams[wdkey] == "":
            self.jobParams[wdkey] = os.getcwd()

        # get the submit options
        self.submitOpts = self.clusterCfg["submit_opts"]

        # get the run cmd options
        self.runCmdOpts = self.clusterCfg["run_cmd"]

        # get the sync job options
        self.syncOpts = self.clusterCfg["sync_job"]

        # analysis_pipeline path in the docker image
        self.pipelinePath = self.clusterCfg["pipeline_path"]

        # create the batch client
        self.batchC = boto3.client('batch',region_name=self.clusterCfg["aws_region"])

    def createDependsList(self, holdid):
        # holdid is either a single job id or a list of job ids (1 or more)
        if  not isinstance(holdid, list):
            holdid = [ holdid ]
        # create a list of jobids that the job depends on
        depends_list = [ {'jobId': id} for id in holdid ]

        return depends_list

    def submitSyncJobs(self, depends_list, jName):
        # submit sync jobs for each each set of 20 jobIDs and return jobids for
        # each sync job
        subIDs = []
        self.printVerbose("2>> submitSyncJobs: job/depends_list: " + jName + "/" + str(depends_list))
        # set the synjobparams
        self.syncOpts["parameters"]["jids"] = str([ d['jobId'] for d in depends_list ])
        # submit sync job in batches of 20
        maxDepends = 20
        noDepends = len(depends_list)
        noSyncJobs = int(math.ceil(noDepends/(maxDepends+1))) + 1
        noDependsLast = noDepends % maxDepends
        if noDependsLast == 0:
            noDependsLast = maxDepends
        if noSyncJobs > maxDepends:
            sys.exit("Error: Too many depend jobs (" + str(noDepends) + ").  Max number is " + str(maxDepends))
        self.printVerbose("2>>>> submitSyncJobs: No. holds/sync jobs/noLast: " + str(noDepends) + "/" + str(noSyncJobs) +
                          "/" + str(noDependsLast))

        for sj in range(noSyncJobs):
            sIndex = sj*maxDepends
            lIndex = sIndex+maxDepends
            if sj == noSyncJobs - 1:
                lIndex = sIndex+noDependsLast
            jobName = 'sync_' + jName + '_' + str(sj)
            self.printVerbose("2>>>> submitSyncJobs: Sumbitting sync job: " + jobName +
                              " depend list[    " + str(sIndex) + "," + str(lIndex) + "] \n\t\t" + str(depends_list[sIndex:lIndex]))
            subid = self.batchC.submit_job(
               jobName = jobName,
               jobQueue = self.syncOpts["submit_opts"]["queue"],
               jobDefinition = self.syncOpts["submit_opts"]["jobdef"],
               parameters = self.syncOpts["parameters"],
               dependsOn = depends_list[sIndex:lIndex])
            subIDs.append(subid)
        # get the job ids for the above sync jobs
        syncDepends_list = [{'jobId': d['jobId']} for d in subIDs]
        # now sumbit a master sync job that waits for the other sync jobs
        # to complete;  just return one jobid that the caller can use
        masterParams = {'msg': 'master job completed', 'jids': jName}
        jobName = 'Mastersync_' + jName
        self.printVerbose("2>> submitSyncJobs: Sumbitting Mastersync job: " + jobName)

        subid = self.batchC.submit_job(
           jobName = jobName,
           jobQueue = self.syncOpts["submit_opts"]["queue"],
           jobDefinition = self.syncOpts["submit_opts"]["jobdef"],
           parameters = masterParams,
           dependsOn = syncDepends_list)
        masterDepend_list = [ {'jobId': subid['jobId']} ]
        self.printVerbose("2>> submitSyncJobs: returning Mastersync return job id list: " + str())

        return masterDepend_list

    def runCmd(self, job_name, cmd, logfile=None):
        runCmdOpts = deepcopy(self.runCmdOpts)
        # set rdriver
        key = "--rdriver"
        if key not in runCmdOpts["params"].keys() or runCmdOpts["params"][key] == "":
            # if use pipelinePath from cfg if defined
            pipelinePath = os.path.dirname(os.path.abspath(cmd[0]))
            baseName = os.path.basename(cmd[0])
            if self.pipelinePath != "":
                pipelinePath = self.pipelinePath
            runCmdOpts["params"][key] = os.path.join(pipelinePath, baseName)

        # set rargs
        key = "--rargs"
        if key not in runCmdOpts["params"].keys() or runCmdOpts["params"][key] == "":
            cs = " ".join(cmd[1:])
            qcs = '"' + cs + '"'
            runCmdOpts["params"][key] = qcs

        # working directory
        key = "--workdir"
        if key not in runCmdOpts["params"].keys() or runCmdOpts["params"][key] == "":
            runCmdOpts["params"][key] = self.jobParams['wd']

        # convert params dict to key-value pair
        sCmdArgs = " ".join([" ".join([key, str(val)]) for key, val in runCmdOpts["params"].items()])

        # create the cmd to spawn
        sCmd = " ".join([runCmdOpts["cmd"], sCmdArgs])

        # redirect stdout/stderr
        if logfile != None:
            sout = sys.stdout
            serr = sys.stderr
            flog = open ( logfile, 'w' )
            sys.stderr = sys.stdout = flog
        # spawn
        print("runCmd: executing " + sCmd)
        process = subprocess.Popen(sCmd, stdout=sys.stdout, stderr=sys.stderr, shell=True)
        status = process.wait()

        # redirect stdout back
        if logfile != "":
            sys.stdout = sout
            sys.stderr = serr
        if status:
            print( "Error running job: " + job_name + " (" + str(status) + ") for command:" )
            print( "\t> " + str(sCmd) )
            sys.exit(2)


    def submitJob(self, job_name, cmd, args=None, holdid=None, array_range=None, request_cores=None, print_only=False, **kwargs):
        self.printVerbose("1>> submitJob: " + job_name + " beginning ...")
        jobParams = deepcopy(self.jobParams)
        submitOpts = deepcopy(self.submitOpts)
        pipelinePath = self.pipelinePath
        # set the R driver and arguments (e.g., -s rcode cfg --chr cn)
        key = "rd"
        # if pipeline path not defined in cfg, then use path from cmd
        if pipelinePath == "":
            pipelinePath = os.path.dirname(os.path.abspath(cmd))
        baseName = os.path.basename(cmd)
        jobParams[key] = os.path.join(pipelinePath, baseName)

        key = "ra"
        if args is None:
            args = []
        jobParams[key] = " ".join(args)

        # assign the log file
        key = "lf"
        t = str(int(time.time()))
        lfile = job_name + "_" + t + ".log"
        jobParams[key] = lfile

        # check for number of cores (1 core = 2 vcpus)
        key = "vcpus"
        if request_cores is not None:
            submitOpts[key] = 2*request_cores

        # get memory limit option
        key = "memory_limits"
        if key in self.clusterCfg.keys():
            # get the memory limits
            memlim = super(AWS_Batch, self).memoryLimit(self.clusterCfg[key], job_name)
            if memlim != None:
                submitOpts["memory"] = memlim

        # holdid is a list of a list of jid dictionaries ([{'jobId': 'jobid string'}, ... ])
        if holdid is not None:
            submitOpts["dependsOn"] = self.createDependsList(holdid)

        # if we're doing a job array equivalent submit multiple jobs; else just submit one job
        submitJobs = []
        if not print_only:
            # with the current limit of jobs that a submitted job can depend upon (20),
            # we'll use of sync jobs to sync if we have more than 20 hold jobs
            if len(self.submitOpts["dependsOn"]) > 0:
                depends_list = self.submitSyncJobs(submitOpts["dependsOn"], job_name)
            else:
                depends_list = submitOpts["dependsOn"]

            if array_range is not None:
                lf = jobParams['lf']
                air = [ int(i) for i in array_range.split( '-' ) ]
                taskList = range( air[0], air[1]+1 )
                self.printVerbose("1>>>> submitJob: " + job_name + " is an array job - no. tasks: " + str(len(taskList)))
                self.printVerbose("1>>>> submitJob: Submitting array jobs with depend list: " + str(depends_list))

                for t in taskList:
                    # add an environmnent for equiv to taskid in sge
                    submitOpts["env"] = [ { "name": "SGE_TASK_ID",
                                            "value": str(t) } ]
                    jobParams['lf'] = lf + '.' + str(t)
                    # create a jobname based on the job_name submitted and the task id
                    subOut = self.batchC.submit_job(
                       jobName = job_name + "_" + str(t),
                       jobQueue = submitOpts["queue"],
                       jobDefinition = submitOpts["jobdef"],
                       parameters = jobParams,
                       dependsOn = depends_list,
                       containerOverrides = {
                          "vcpus": submitOpts["vcpus"],
                          "memory": submitOpts["memory"],
                          "environment": submitOpts["env"]
                       }
                    )
                    # return from batch submit_job is a dict of the form: { 'jobId': 'xxx-yyy', 'jobName': 'myjobName' },
                    # which we'll return as a list of submit jobs dict
                    submitJobs.append( subOut )
                # conver to depend list
                arrayDepends_list = [ {'jobId': d['jobId']} for d in submitJobs]
                # sumbit a sync job for the array
                msJobName = job_name+'_tasks_'+array_range
                self.printVerbose("1>>>> submitJob: Submitting Mastersync job " + msJobName +
                                  " for array jobs")
                depends_list = self.submitSyncJobs(arrayDepends_list, job_name+'_tasks'+array_range)
                jobid = depends_list[0]['jobId']
            else:
                self.printVerbose("1>>>> submitJob: Submitting single job " + job_name + " with depend list: " + str(depends_list))

                subOut = self.batchC.submit_job(
                   jobName = job_name,
                   jobQueue = submitOpts["queue"],
                   jobDefinition = submitOpts["jobdef"],
                   parameters = jobParams,
                   dependsOn = depends_list,
                   containerOverrides = {
                      "vcpus": submitOpts["vcpus"],
                      "memory": submitOpts["memory"],
                      "environment": submitOpts["env"]
                   }
                )
                # return from batch submit_job jus the jobId
                jobid = subOut['jobId']
        else:
            batchCmd = "\n\tjobName = " + job_name
            batchCmd = batchCmd + ", jobQueue = " + submitOpts["queue"]
            batchCmd = batchCmd + ", jobDef = " + submitOpts["jobdef"]
            batchCmd = batchCmd + ", memory = " + str(submitOpts["memory"])
            batchCmd = batchCmd + ", vcpus = " + str(submitOpts["vcpus"])
            batchCmd = batchCmd + "\n\tparams = " + str(jobParams)
            if array_range is not None:
                ct = "array job " + array_range
            else:
                ct = "single job"
            ht = str(submitOpts["dependsOn"])
            jobid = "111-222-333-print_only-" +  job_name

            print( "\njob name: " + job_name + " job type: " + ct + "\tdepends: " + ht + "\tjobid:" + jobid +
                   "\nbatch.submit_job : " + batchCmd )
            # jids from batch is a dict of the form: { 'jobId': 'xxx-yyy', 'jobName': 'myjobName' }

        # return the job id (either from the single job or array job)
        self.printVerbose("1>> submitJob: " + job_name + " returning jobid: " + jobid)

        return jobid


class SGE_Cluster(Cluster):

    def __init__(self, cluster_file=None, verbose=True):
        self.class_name = self.__class__.__name__
        super(SGE_Cluster, self).__init__(cluster_file, verbose)

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
        # get memory limit option
        key = "memory_limits"
        if key in self.clusterCfg.keys():
            memlim = super(SGE_Cluster, self).memoryLimit(self.clusterCfg[key], kwargs["job_name"])
            if memlim != None:
                submit_opts["-l"] = "h_vmem="+str(memlim)+"G"

        jobid = self.executeJobCmd(subOpts, **kwargs)
        return jobid

    def executeJobCmd(self, submitOpts, **kwargs):
        # update sge opts
        submitOpts["-N"] = kwargs["job_name"]

        key = "holdid"
        if key in kwargs and kwargs["holdid"] != []:
            if isinstance(kwargs["holdid"], str):
                kwargs["holdid"] = [kwargs["holdid"]]
            submitOpts["-hold_jid"] =  ",".join(kwargs["holdid"])

        key = "array_range"
        if key in kwargs:
            submitOpts["-t"] = kwargs["array_range"]

        key = "request_cores"
        if key in kwargs and kwargs[key] != None and kwargs[key] != "1":
            submitOpts["-pe"] = self.clusterCfg["parallel_env"] + " " + kwargs["request_cores"]

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
        process = subprocess.Popen(sub_cmd, shell=True, stdout=subprocess.PIPE)
        pipe = process.stdout
        sub_out = pipe.readline()
        jobid = sub_out.strip(' \t\n\r')

        if "array_range" in kwargs:
            jobid = jobid.split(".")[0]
        print("Submitting job " + jobid + " (" + kwargs["job_name"] + ")")

        return jobid


class UW_Cluster(SGE_Cluster):

    def __init__(self, cluster_file=None, verbose=False):
        self.class_name = self.__class__.__name__
        super(UW_Cluster, self).__init__(cluster_file, verbose)


class AWS_Cluster(SGE_Cluster):

    def __init__(self, cluster_file=None, verbose=False):
        self.class_name = self.__class__.__name__
        super(AWS_Cluster, self).__init__(cluster_file, verbose)

    def submitJob(self, **kwargs):
        # currently, no email on aws
        kwargs["email"] = None
        jobid = super(AWS_Cluster, self).submitJob(**kwargs)
        return jobid


class ClusterFactory(object):

    @staticmethod
    def createCluster(cluster_type, cluster_file, verbose):
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
