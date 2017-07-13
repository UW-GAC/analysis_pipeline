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
    def __init__(self, cluster_file=None, verbose=False):
        self.verbose = verbose
        self.class_name = self.__class__.__name__
        self.openClusterCfg("./cluster_cfg.json", cluster_file, verbose)

    def openClusterCfg(self, stdCfgFile, optCfgFile, verbose):
        # get the standard cluster cfg
        pipePath = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.clusterfile =  os.path.join(pipePath, "./cluster_cfg.json")
        self.printVerbose("0>> openClusterCfg: Reading internal cfg file: " + self.clusterfile)

        with open(self.clusterfile) as cfgFileHandle:
            clusterCfg= json.load(cfgFileHandle)
        key = "debug"
        debugCfg = False
        if key in clusterCfg:
            if clusterCfg[key] == 1:
                debugCfg = True
        self.clusterCfg = clusterCfg["cluster_types"][self.class_name]
        if debugCfg:
            print("0>>> Dump of standard cluster cfg ... \n")
            print json.dumps(self.clusterCfg, indent=3, sort_keys=True)
        if optCfgFile != None:
            self.printVerbose("0>> openClusterCfg: Reading user cfg file: " + self.clusterfile)

            with open(optCfgFile) as cfgFileHandle:
                clusterCfg = json.load(cfgFileHandle)
            optCfg = clusterCfg["cluster_types"][self.class_name]
            if debugCfg:
                print("0>>> Dump of optional cluster cfg ... \n")
                print json.dumps(optCfg, indent=3, sort_keys=True)
            # update
            self.clusterCfg = update(self.clusterCfg, optCfg)
            if debugCfg:
                print("0>>> Dump of updated cluster cfg ... \n")
                print json.dumps(self.clusterCfg, indent=3, sort_keys=True)

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
        # get the queue
        self.queue = self.clusterCfg["queue"]

        # create the batch client
        self.batchC = boto3.client('batch',region_name=self.clusterCfg["aws_region"])

    def getIDsAndNames(self, submitHolds):
        # for the submit holds, return a dictionary of all job names in a single string
        # and a list of all job ids
        nlist = [name for d in submitHolds for name in d]
        jobnames = "_".join(nlist)
        jobids = [id for d in submitHolds for il in d.values() for id in il]
        return {'jobnames': jobnames, 'jobids': jobids}

    def submitSyncJobs(self, job_name, submitHolds):
        # create a list of {'jobId': jobid} compatible with batch submit job associated with the
        # submit holds. if no. of jobids > 20, create two or more sync jobs and return those jobids
        holds = self.getIDsAndNames(submitHolds)
        jids = holds['jobids']
        hold_jnames = holds['jobnames']
        dependsList = [{'jobId': jid} for jid in jids]
        self.printVerbose("\t2> submitSyncJobs: job " + job_name + " depends on " + hold_jnames + " with " + str(len(jids)) + " job ids")
        maxDepends = 20
        if len(jids)> maxDepends:
            self.printVerbose("\t2> submitSyncJobs: job " + job_name + " - creating intemediary sync jobs ...")
            # set the synjobparams
            self.syncOpts["parameters"]["jids"] = str(jids)
            # submit sync job in batches of 20
            maxDepends = 20
            noDepends = len(jids)
            noSyncJobs = int(math.ceil(noDepends/(maxDepends+1))) + 1
            noDependsLast = noDepends % maxDepends
            if noDependsLast == 0:
                noDependsLast = maxDepends
            if noSyncJobs > maxDepends:
                sys.exit("Error: Too many hold jobs to sync_ (" + str(noDepends) + ").  Max number of sync jobs is " + str(maxDepends))
            self.printVerbose("\t\t2>> submitSyncJobs: No. holds/sync jobs/noLast: " + str(noDepends) + "/" + str(noSyncJobs) +
                              "/" + str(noDependsLast))
            syncDepends_list = []
            for sj in range(noSyncJobs):
                sIndex = sj*maxDepends
                lIndex = sIndex+maxDepends
                if sj == noSyncJobs - 1:
                    lIndex = sIndex+noDependsLast
                jobName = job_name + '_DependsOn_' + hold_jnames + '_' + str(sj)
                self.printVerbose("\t\t2>> submitSyncJobs: Sumbitting sync job: " + jobName +
                                  " depend list[    " + str(sIndex) + "," + str(lIndex) + "] \n\t\t\t" + str(dependsList[sIndex:lIndex]))
                subid = self.batchC.submit_job(
                   jobName = jobName,
                   jobQueue = self.queue,
                   jobDefinition = self.syncOpts["submit_opts"]["jobdef"],
                   parameters = self.syncOpts["parameters"],
                   dependsOn = dependsList[sIndex:lIndex])
                syncDepends_list.append({'jobId': subid['jobId']})
            dependsList = syncDepends_list

        self.printVerbose("\t2> submitSyncJobs: job " + job_name + " will depend on the job ids:\n\t\t" + str(dependsList))
        return dependsList

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
        self.printVerbose("1===== submitJob: job " + job_name + " beginning ...")
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
            args = ["NoArgs"]
        jobParams[key] = " ".join(args)

        # using time set a job id (which is for tracking; not the batch job id)
        trackID = job_name + "_" + str(int(time.time()*100))
        logExt = ".log"

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

        # holdid is a previous submit_id dict {job_name: [job_Ids]}
        if holdid is not None:
            submitHolds = holdid
        else:
            submitHolds = []

        # if we're doing a job array equivalent submit multiple jobs; else just submit one job
        if not print_only:
            # see if there are any holdids
            if len(submitHolds) > 0:
                # process hold ids and return a "dependsOn" list
                submitOpts["dependsOn"] = self.submitSyncJobs(job_name, submitHolds)
            # set the log file name that's common to both single and array jobs
            if len(submitOpts["dependsOn"]) > 0:
                self.printVerbose("\t1> submitJob: " + job_name + " depends on " + self.getIDsAndNames(submitHolds)['jobnames'])
            else:
                self.printVerbose("\t1> submitJob: " + job_name + " does not depend on other jobs" )

            log_prefix = trackID
            if array_range is not None:
                air = [ int(i) for i in array_range.split( '-' ) ]
                taskList = range( air[0], air[1]+1 )
                self.printVerbose("\t1> submitJob: " + job_name + " is an array job - no. tasks: " + str(len(taskList)))
                hold_name = job_name + "_" + trackID
                job_ids=[]
                for t in taskList:
                    # add an environmnent for equiv to taskid in sge; and jobid
                    submitOpts["env"].append( { "name": "SGE_TASK_ID",
                                                "value": str(t) } )
                    submitOpts["env"].append( { "name": "JOB_ID",
                                                "value": trackID } )
                    jobParams['lf'] = log_prefix + '.task_' + str(t) + logExt
                    # create a jobname based on the job_name submitted and the task id
                    subOut = self.batchC.submit_job(
                       jobName = job_name + "_" + str(t),
                       jobQueue = self.queue,
                       jobDefinition = submitOpts["jobdef"],
                       parameters = jobParams,
                       dependsOn = submitOpts["dependsOn"],
                       containerOverrides = {
                          "vcpus": submitOpts["vcpus"],
                          "memory": submitOpts["memory"],
                          "environment": submitOpts["env"]
                       }
                    )
                    # append to the list of job ids
                    job_ids.append(subOut['jobId'])
                # create the submit id dictionary
                submit_id = {job_name: job_ids}
            else:
                self.printVerbose("\t1> submitJob: " + job_name + " is a single job")
                jobParams['lf'] = log_prefix + logExt
                submitOpts["env"].append( { "name": "JOB_ID",
                                            "value": trackID } )
                subOut = self.batchC.submit_job(
                   jobName = job_name,
                   jobQueue = self.queue,
                   jobDefinition = submitOpts["jobdef"],
                   parameters = jobParams,
                   dependsOn = submitOpts["dependsOn"],
                   containerOverrides = {
                      "vcpus": submitOpts["vcpus"],
                      "memory": submitOpts["memory"],
                      "environment": submitOpts["env"]
                   }
                )
                # return the "submit_id" which is a list of dictionaries
                submit_id = {job_name: [subOut['jobId']]}
        else:
            log_prefix = trackID
            jobParams['lf'] = log_prefix + "_printonly" + logExt
            print("Submit job: " + job_name)
            if array_range is not None:
                print("\tjob is an array job: " + array_range)
            else:
                print("\tjob is a single job")
            print("\tbatch queue: " + self.queue)
            print("\tjob definition: " + submitOpts["jobdef"])
            print("\tjob memory: " + str(submitOpts["memory"]))
            print("\tjob vcpus: " + str(submitOpts["vcpus"]))
            print("\tjob params: \n\t\t" + str(jobParams))
            jobid = "111-222-333-print_only-" +  job_name
            submit_id = {job_name: [jobid]}
            print("\tsubmit_id: " + str(submit_id))

        # return the job id (either from the single job or array job)
        self.printVerbose("\t1> submitJob: " + job_name + " returning submit_id: " + str(submit_id))

        return submit_id


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
                subOpts["-l"] = "h_vmem="+str(memlim)+"G"

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
