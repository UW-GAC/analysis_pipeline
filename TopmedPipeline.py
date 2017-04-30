"""Utility functions for TOPMed pipeline"""

import sys
import csv
import subprocess
from copy import deepcopy
import getpass
import time
import os
import json

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

# parent class to represent a compute cluster environment
class Cluster(object):
    """ """
    # constructor
    def __init__(self, submit_cmd, cluster_file=None):
        self.class_name = self.__class__.__name__
        # get the standard cluster cfg
        pipePath = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.clusterfile =  os.path.join(pipePath, "./cluster_cfg.json")
        with open(self.clusterfile) as cfgFileHandle:
            clusterCfg= json.load(cfgFileHandle)
        self.clusterCfg = clusterCfg["cluster_types"][self.class_name]
        if cluster_file != None:
            with open(cluster_file) as cfgFileHandle:
                clusterCfg = json.load(cfgFileHandle)
            optCfg = clusterCfg["cluster_types"][self.class_name]
            # update
            self.clusterCfg.update(optCfg)

    def getClusterCfg(self):
        return self.clusterCfg

    def memoryLimit(self, memLimits, job_name):
        memlim = None
        a =[ akey for akey in memLimits.keys() if job_name.startswith( akey ) ]
        if len(a):
            # just find the first match to job_name
            memlim = memLimits[a[0]]
        return memlim

class AWS_Batch(Cluster):

    def __init__(self, cluster_file=None):
        self.class_name = self.__class__.__name__
        super(AWS_Batch, self).__init__(cluster_file)
        self.clusterCfg = super(AWS_Batch, self).getClusterCfg()

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

    def getSyncJobs(self, fullJobIDs, jName):
        # submit sync jobs for each each set of 20 jobIDs and return jobids for
        # each sync job
        holdJobNames = [fullJobIDs[i]['jobName'] for i in range(len(fullJobIDs))]
        holdJobIDs = [fullJobIDs[i]['jobID'] for i in range(len(fullJobIDs))]
        syncIDs = []
        # set the synjobparams
        self.syncOpts["parameters"]["hjn"] = holdJobNames
        # submit sync job in batches of 20
        maxHolds = 20
        noHolds = len(fullJobIDs)
        noSyncJobs = math.ceil(noHolds/maxHolds)
        noHoldsLast = noHolds % maxHolds
        if noHoldsLast == 0:
            noHoldsLast = maxHolds
        if noSyncJobs > maxHolds:
            sys.exit("Error: Too many sync jobs (" + str(numberSyncJobs) + ")")
        for sj in range(noSyncJobs):
            if sj == noSyncJobs - 1:
                noHolds = noHoldsLast
            else:
                noHolds = maxHolds
            holdIndices = range(sj*maxHolds,sj*maxHolds+noholds)

            noHoldJobs = len(sj)
            syncID = self.batchC.submit_job(
               jobName = 'sync_' + jName + '_' + str(sj),
               jobQueue = self.syncOpts["submit_opts"]["queue"],
               jobDefinition = self.syncOpts["submit_opts"]["jobdef"],
               parameters = syncJobParams,
               dependsOn = [holdJobNames[n] for n in holdIndices ])
            syncIDs.append(syncID)

        return syncIDs

    def runCmd(self, job_name, cmd, logfile=None):
        # set rdriver
        key = "--rdriver"
        if key not in self.runCmdOpts["params"].keys() or self.runCmdOpts["params"][key] == "":
            # if use pipelinePath from cfg if defined
            pipelinePath = os.path.dirname(os.path.abspath(cmd[0]))
            baseName = os.path.basename(cmd[0])
            if self.pipelinePath != "":
                pipelinePath = self.pipelinePath
            self.runCmdOpts["params"][key] = os.path.join(pipelinePath, baseName)

        # set rargs
        key = "--rargs"
        if key not in self.runCmdOpts["params"].keys() or self.runCmdOpts["params"][key] == "":
            cs = " ".join(cmd[1:])
            qcs = '"' + cs + '"'
            self.runCmdOpts["params"][key] = qcs

        # working directory
        key = "--workdir"
        if key not in self.runCmdOpts["params"].keys() or self.runCmdOpts["params"][key] == "":
            self.runCmdOpts["params"][key] = self.jobParams['wd']

        # convert params dict to key-value pair
        sCmdArgs = " ".join([" ".join([key, str(val)]) for key, val in self.runCmdOpts["params"].items()])

        # create the cmd to spawn
        sCmd = " ".join([self.runCmdOpts["cmd"], sCmdArgs])

        # redirect stdout/stderr
        if logfile != None:
            sout = sys.stdout
            serr = sys.stderr
            flog = open ( logfile, 'w' )
            sys.stderr = sys.stdout = flog
        # spawn
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


    def submitJob(self, job_name, cmd, args=None, holdid=None, array_range=None, print_only=False, **kwargs):

        # set the R driver and arguments (e.g., -s rcode cfg --chr cn)
        key = "rd"
        # if use pipelinePath from cfg if defined
        pipelinePath = os.path.dirname(os.path.abspath(cmd))
        baseName = os.path.basename(cmd)
        if self.pipelinePath != "":
            pipelinePath = self.pipelinePath
        self.jobParams[key] = os.path.join(pipelinePath, baseName)

        key = "ra"
        if args is None:
            args = []
        self.jobParams[key] = " ".join(args)

        # assign the log file
        key = "lf"
        if key not in self.jobParams or self.jobParams[key] == "":
            t = str(int(time.time()))
            lfile = job_name + "_" + t + ".log"
            self.jobParams[key] = lfile

        # get memory limit option
        key = "memory_limits"
        if key in self.clusterCfg.keys():
            # get the memory limits
            memlim = super(AWS_Batch, self).memoryLimit(self.clusterCfg[key], job_name)
            if memlim != None:
                self.submitOpts["memory"] = memlim

        # holdid is a list of a list of jid dictionaries
        if holdid is not None and holdid[0] != []:
            self.submitOpts["dependsOn"] = holdid[0]

        # if we're doing a job array equivalent submit multiple jobs; else just submit one job
        idsonly = []
        jids = []
        if not print_only:
            if len(self.submitOpts["dependsOn"]):
                syncJobIDs = self.getSyncJobs(self.submitOpts["dependsOn"], job_name)
            else:
                syncJobIDs = []
            if array_range is not None:
                lf = self.defParams['lf']
                air = [ int(i) for i in array_range.split( '-' ) ]
                for t in range( air[0], air[1]+1 ):
                    # add an environmnent for equiv to taskid in sge
                    self.defEnv = [ { "name": "SGE_TASK_ID",
                                      "value": str(t) } ]
                    self.defParams['lf'] = lf + '.' + str(t)

                    subOut = self.batchC.submit_job(
                       jobName = job_name,
                       jobQueue = self.submitOpts["queue"],
                       jobDefinition = self.submitOpts["jobdef"],
                       parameters = self.jobParams,
                       dependsOn = syncJobIDs,
                       containerOverrides = {
                          "vcpus": self.submitOpts["vcpus"],
                          "memory": self.submitOpts["memory"],
                          "environment": self.submitOpts["env"]
                       }
                    )
                    # append the jid dictionary to list
                    jids.append( subOut )
                    idsonly.append( subOut["jobId"] )
            else:
                subOut = self.batchC.submit_job(
                   jobName = job_name,
                   jobQueue = self.submitOpts["queue"],
                   jobDefinition = self.submitOpts["jobdef"],
                   parameters = self.jobParams,
                   dependsOn = syncJobIDs,
                   containerOverrides = {
                      "vcpus": self.submitOpts["vcpus"],
                      "memory": self.submitOpts["memory"],
                      "environment": self.submitOpts["env"]
                   }
                )
                # append the jid dictionary to list
                jids.append( subOut )
                idsonly.append( subOut["jobId"] )
            print( "Job Name: " + job_name + " Job IDs: " + str(idsonly) )
        else:
            batchCmd = "\n\tjobName = " + job_name
            batchCmd = batchCmd + ", jobQueue = " + self.submitOpts["queue"]
            batchCmd = batchCmd + ", jobDef = " + self.submitOpts["jobdef"]
            batchCmd = batchCmd + ", memory = " + str(self.submitOpts["memory"])
            batchCmd = batchCmd + ", vcpus = " + str(self.submitOpts["vcpus"])
            batchCmd = batchCmd + "\n\tparams = " + str(self.jobParams)
            if array_range is not None:
                ct = "array job " + array_range
            else:
                ct = "single job"
            if holdid is None or holdid[0] == []:
                ht = "no hold"
            else:
                ht = holdid[0]

            print( "\njob name: " + job_name + " job type: " + ct + "\tht: " + str(ht) + "\nbatch.submit_job : " + batchCmd)
            jids = [ { 'jobId': 'print_only-1/' + job_name },
                     { 'jobId': 'print_only-2/' + job_name } ]

            # return a list of jid dictionaries [ {'jobId': string_of_jobid}, ...]
        return jids

class SGE_Cluster(Cluster):

    def __init__(self, cluster_file=None):
        self.class_name = self.__class__.__name__
        super(SGE_Cluster, self).__init__(cluster_file)
        self.clusterCfg = super(UW_Cluster, self).getClusterCfg()

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
        # set the job cmd
        kwargs["job_cmd"] = self.clusterCfg["submit_cmd"]
        # get memory limit option
        key = "memory_limits"
        if key in self.clusterCfg.keys():
            memlim = super(SGE_Cluster, self).memoryLimit(self.clusterCfg[key], kwargs["job_name"])
            if memlim != None:
                self.submit_opts["-l"] = "h_vmem="+str(memlim)+"G"

        jobid = self.executeJobCmd(self.clusterCfg["submit_opts"], **kwargs)
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
        if key in kwargs:
            submitOpts["-pe"] = "local " + kwargs["request_cores"]

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

        process = subprocess.Popen(sub_cmd, shell=True, stdout=subprocess.PIPE)
        pipe = process.stdout
        sub_out = pipe.readline()
        jobid = sub_out.split()[2]

        if verbose:
            print sub_out

        key = "array_range"
        if key in kwargs:
            jobid = jobid.split(".")[0]

        return jobid

class UW_Cluster(SGE_Cluster):

    def __init__(self, cluster_file=None):
        self.class_name = self.__class__.__name__
        super(UW_Cluster, self).__init__(cluster_file)
        self.clusterCfg = super(UW_Cluster, self).getClusterCfg()


class AWS_Cluster(SGE_Cluster):

    def __init__(self, cluster_file=None):
        self.class_name = self.__class__.__name__
        super(AWS_Cluster, self).__init__(cluster_file)
        self.clusterCfg = super(AWS_Cluster, self).getClusterCfg()

    def submitJob(self, **kwargs):
        # currently, no email on aws
        kwargs["email"] = None
        jobid = super(AWS_Cluster, self).submitJob(**kwargs)
        return jobid

class ClusterFactory(object):

    @staticmethod
    def createCluster(cluster_type, cluster_file):
        allSubClasses = getAllSubclasses(Cluster)
        for subclass in allSubClasses:
            if subclass.__name__ == cluster_type:
                return subclass(cluster_file)
        raise Exception("unknown cluster type: " + cluster_type + "!")

def getAllSubclasses(base):
    all_subclasses = []
    for subclass in base.__subclasses__():
        all_subclasses.append(subclass)
        all_subclasses.extend(getAllSubclasses(subclass))
    return all_subclasses
