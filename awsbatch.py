#! /usr/bin/env python3
# submit a single analysis from the analysis pipeline (e.g., ld_pruning)
from __future__ import division
import      os
import      getpass
import      time
import      json
import      math
import      collections
import      sys
from        copy   import deepcopy
from        argparse import ArgumentParser
from        datetime import datetime, timedelta

try:
    import boto3
except ImportError:
    batchSupport = False
else:
    batchSupport = True

# init globals
fileversion = '1.0'
msgErrPrefix = '>>> Error: '
msgInfoPrefix = '>>> Info: '
debugPrefix = '>>> Debug: '
verbose = False

batchClient = None
defApath = "/usr/local/analysis_pipeline"
defBaseCmd = "runRscript.sh"
subParams = {   'profile': None,
                'cluster_file': None,
                'clustercfg': None,
                'jobname': None,
                'cmd': None,
                'args': None,
                'holdid': None,
                'array_range': None,
                'request_cores': None,
                'print_only': None,
                'verbose': None,
                'maxmem': None,
                'workdir': None,
                'queue': None,
                'jobdef': None,
                'apath': None,
                'infofile': None,
                'analysislog': None
            }

def getBatchClient(profile_a):
    global batchClient

    if not batchSupport:
        pError('getBatchClient: AWS Batch is not supported (boto3 is not available)')
        sys.exit(2)
    if batchClient == None:
        # create the batch client
        try:
            session = boto3.Session(profile_name = profile_a)
            batchClient = session.client('batch')
        except Exception as e:
            pError('boto3 session or client exception ' + str(e))
            sys.exit(2)
    return batchClient

# cluster configuration is read from json into nested dictionaries
# regular dictionary update loses default values below the first level
# https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
def updatecfg(d, u):
    ld = deepcopy(d)
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            if len(v) == 0:
                ld[k] = u[k]
            else:
                r = updatecfg(d.get(k, {}), v)
                ld[k] = r
        else:
            ld[k] = u[k]
    return ld

def memoryLimit(job_name, clusterCfg):
    memlim = None
    memLimits = clusterCfg["memory_limits"]
    if memLimits is None:
        return memlim
    jobMem = [ v for k,v in memLimits.items() if job_name.find(k) != -1 ]
    if len(jobMem):
        # just find the first match to job_name
        memlim = jobMem[0]
    return memlim

def getClusterCfg(a_stdcfg, a_optcfg, a_cfgversion):
    # get the standard cluster cfg
    with open(a_stdcfg) as cfgFileHandle:
        clustercfg = json.load(cfgFileHandle)
    # check version
    key = "version"
    if key in clustercfg:
        if clustercfg[key] != a_cfgversion:
            pError( "Incorrect version of " + a_stdcfg + ": " + clustercfg[key] +
                   "(" + a_cfgversion + ")")
            sys.exit(2)
    else:
        pError( "Version key missing in " + a_stdcfg )
        sys.exit(2)
    key = "debug"
    debugCfg = False
    if key in clustercfg:
        if clustercfg[key] == 1:
            debugCfg = True
    key = "configuration"
    if key in clustercfg:
        clusterconfig = clustercfg[key]
        if debugCfg:
            pInfo("Dump of " + clustercfg["name"] + " ... \n")
            print(json.dumps(clusterconfig, indent=3, sort_keys=True))
        if a_optcfg != None:
            pDebug("Option cluster cfg file: " + a_optcfg)

            with open(a_optcfg) as cfgFileHandle:
                optcfg = json.load(cfgFileHandle)
            optconfiguration = optcfg["configuration"]
            if debugCfg:
                pDebug("Dump of " + optcfg["name"] + " ... \n")
                print(json.dumps(optconfiguration, indent=3, sort_keys=True))
            # update
            clusterconfig = updatecfg(clusterconfig, optconfiguration)
            if debugCfg:
                pDebug("Dump of updated cluster cfg ... \n")
                print(json.dumps(clusterconfig, indent=3, sort_keys=True))
    return clusterconfig

def getIDsAndNames(submitHolds):
    # for the submit holds, return a dictionary of all job names in a single string
    # and a list of all job ids
    nlist = [name for d in submitHolds for name in d]
    maxLen = 1
    if len(nlist) > maxLen:
        nlist = nlist[:maxLen]
    jobnames = "_".join(nlist) + "_more"
    jobids = [id for d in submitHolds for il in list(d.values()) for id in il]
    return {'jobnames': jobnames, 'jobids': jobids}

def submitSyncJobs(job_name, submitHolds, clustercfg, queue):
    # create a list of {'jobId': jobid} compatible with batch submit job associated with the
    # submit holds. if no. of jobids > 20, create two or more sync jobs and return those jobids
    holds = getIDsAndNames(submitHolds)
    jids = holds['jobids']
    hold_jnames = holds['jobnames']
    dependsList = [{'jobId': jid} for jid in jids]
    syncOpts = clustercfg["sync_job"]
    pDebug("\t2> submitSyncJobs: job " + job_name + " depends on " + hold_jnames + " with " + str(len(jids)) + " job ids")
    maxDepends = 20
    if len(jids)> maxDepends:
        pDebug("\t2> submitSyncJobs: job " + job_name + " - creating intemediary sync jobs ...")
        # set the synjobparams
        syncOpts["parameters"]["jids"] = str(jids)
        # submit sync job in batches of 20
        maxDepends = 20
        noDepends = len(jids)
        noSyncJobs = int(math.ceil((noDepends/(maxDepends+1)))) + 1
        noDependsLast = noDepends % maxDepends
        if noDependsLast == 0:
            noDependsLast = maxDepends
        if noSyncJobs > maxDepends:
            sys.exit("Error: Too many hold jobs to sync_ (" + str(noDepends) + ").  Max number of sync jobs is " + str(maxDepends))
        pDebug("\t\t2>> submitSyncJobs: No. holds/sync jobs/noLast: " + str(noDepends) + "/" + str(noSyncJobs) +
                          "/" + str(noDependsLast))
        # get the batch client
        batchC = getBatchClient(clustercfg["aws_profile"])

        syncDepends_list = []
        for sj in range(noSyncJobs):
            sIndex = sj*maxDepends
            lIndex = sIndex+maxDepends
            if sj == noSyncJobs - 1:
                lIndex = sIndex+noDependsLast
            jobName = job_name + '_DependsOn_' + hold_jnames + '_' + str(sj)
            pDebug("\t\t2>> submitSyncJobs: Sumbitting sync job: " + jobName +
                              " depend list[    " + str(sIndex) + "," + str(lIndex) + "] \n\t\t\t" + str(dependsList[sIndex:lIndex]))
            subid = batchC.submit_job(
               jobName = jobName,
               jobQueue = queue,
               jobDefinition = syncOpts["submit_opts"]["jobdef"],
               parameters = syncOpts["parameters"],
               dependsOn = dependsList[sIndex:lIndex])
            syncDepends_list.append({'jobId': subid['jobId']})
        dependsList = syncDepends_list

    pDebug("\t2> submitSyncJobs: job " + job_name + " will depend on the job ids:\n\t\t" + str(dependsList))
    return dependsList

def submitjob(a_submitParams):
    # get all the parameters and update clustercfg
    # job name
    job_name = a_submitParams["jobname"]
    if job_name == None:
        pError("job_name has not been specified.")
        sys.exit(2)
    lmsg = "Job: " + job_name

    cluster_file = a_submitParams["cluster_file"]    # for print_only
    # get the cluster cfg dict
    clustercfg = deepcopy(a_submitParams["clustercfg"])
    # get the job parameters
    jobParams = clustercfg["job_parameters"]
    # get the submit options
    submitOpts = clustercfg["submit_opts"]
    # get the sync opts
    syncOpts = clustercfg["sync_job"]
    # profile
    profile = a_submitParams["profile"]
    if profile != None:
        clustercfg["aws_profile"] = profile
    # array job
    arrayJob = False
    lmsg_array = "no"
    array_range = a_submitParams["array_range"]
    if array_range is not None:
        air = [ int(i) for i in array_range.split( '-' ) ]
        taskList = list(range( air[0], air[len(air)-1]+1))
        noJobs = len(taskList)
        if noJobs > 1:
            arrayJob = True
            envName = "FIRST_INDEX"
        else:
            envName = "SGE_TASK_ID"
        # set env variable appropriately
        key = "env"
        if key in submitOpts:
            submitOpts["env"].append( { "name": envName,
                                        "value": str(taskList[0]) } )
        else:
            submitOpts["env"] = [ { "name": envName,
                                    "value": str(taskList[0]) } ]
        lmsg_array = str(noJobs)
    lmsg += " / array: " + lmsg_array
    # base cmd or r drive
    base_cmd = a_submitParams["cmd"]
    apath = a_submitParams["apath"]
    if base_cmd == None:
        base_cmd = defBaseCmd
    if apath == None:
        apath = defApath
    full_cmd = apath + "/" + base_cmd
    if not os.path.isfile(full_cmd):
        pError("Full cmd " + full_cmd + " does not exist")
        sys.exit(2)
    jobParams["rd"] = full_cmd
    # cmd args
    cmd_args = a_submitParams["args"]
    if cmd_args == None:
        pError("Cmd args have not been specified.")
        sys.exit(2)
    jobParams["ra"] = cmd_args
    # working directory
    workdir = a_submitParams["workdir"]
    key = "wd"
    if workdir != None:
        jobParams[key] = workdir
    else:
        jobParams[key] = os.getenv('PWD')
    if not os.path.isdir(jobParams[key]):
        pError("Work directory " + jobParams[key] + " does not exist")
        sys.exit(2)
    # no. cores
    request_cores = a_submitParams["request_cores"]
    key = "vcpus"
    lmsg_vcpus = str(submitOpts[key])
    if request_cores != None:
        ncl = request_cores.split("-")
        nci = int(ncl[len(ncl)-1])
        ncs = str(nci)
        submitOpts[key] = nci
        key2 = "env"
        if key2 in submitOpts:
            submitOpts[key2].append( { "name": "NSLOTS",
                                       "value": ncs } )
            submitOpts[key2].append( { "name": "MAX_NUM_THREADS",
                                       "value": ncs } )
        else:
            submitOpts[key2]=[ { "name": "NSLOTS",
                                 "value": ncs } ]
            submitOpts[key2].append( { "name": "MAX_NUM_THREADS",
                                       "value": ncs } )
        lmsg_vcpus = ncs
    lmsg += " / cores: " + lmsg_vcpus
    # set queue
    queue = a_submitParams["queue"]
    if queue == None:
        queue = clustercfg["queue"]
    else:
        clustercfg["queue"] = queue
    # max mem
    maxmem = a_submitParams["maxmem"]               # optional
    key1 = "memory"
    if maxmem != None:
        submitOpts[key1] = int(maxmem)
    else:
        key2 = "memory_limits"
        if key2 in list(clustercfg.keys()):
            memlim = memoryLimit(job_name, clustercfg)
            if memlim != None:
                submitOpts[key1] = memlim
    lmsg += " / mem: " + str(submitOpts[key1])
    # hold ids
    holdid = a_submitParams["holdid"]               # none for interactive
    if holdid is not None:
        submitHolds = holdid
    else:
        submitHolds = []
    # job def
    jobdef = a_submitParams["jobdef"]
    if jobdef != None:
        submitOpts["jobdef"] = jobdef
    # get/set cluster cfg attributes
    retryStrategy = clustercfg["retryStrategy"]
    # using time set a job id (which is for tracking; not the batch job id)
    trackID = job_name + "_" + str(int(time.time()*100))
    # environment variables for job id/track id
    key = "env"
    if key in submitOpts:
        submitOpts[key].append( { "name": "JOB_ID",
                                  "value": trackID } )
    else:
        submitOpts[key]=[ { "name": "JOB_ID",
                            "value": trackID } ]
    # update the log message
    a_submitParams["logmsg"] = lmsg

    print_only = a_submitParams["print_only"]
    global verbose
    verbose = a_submitParams["verbose"]
    # get the batch client
    batchC = getBatchClient(clustercfg["aws_profile"])
    # if print_only, just print out the submit command compatible with interactive submit
    if print_only:
        scmd = apath + "/" + os.path.basename(__name__) + ".py"
        if cluster_file != None:
            scmd += " --cluster_file " + cluster_file
        if workdir != None:
            scmd += " --workdir " + workdir
        scmd += " -j " + job_name
        if jobdef != None:
            scmd += " --jobdef " + jobdef
        if a_submitParams["cmd"] != None:
            scmd += " --basecmd " + base_cmd
        if a_submitParams["apath"] != None:
            scmd += " --apath " + apath

        scmd += " --param " + '"' + cmd_args + '"'
        if array_range != None:
            scmd += " --array " + array_range
        if request_cores != None:
            scmd += " --nocores " + str(request_cores)
        if maxmem != None:
            scmd += " --maxmem " + maxmem
        if queue != None:
            scmd += " -q " + queue
        if profile != None:
            scmd += " --profile " + profile
        pInfo("Job: " + job_name + " - awsbatch submit job command: \n" + scmd + "\n")
        subOut = {'jobName': job_name, 'jobId': "000000"}
    else:
        if len(submitHolds) > 0:
            # process hold ids and return a "dependsOn" list
            submitOpts["dependsOn"] = submitSyncJobs(job_name, submitHolds, clustercfg, queue)
        # set the log file name that's common to both single and array jobs
        if len(submitOpts["dependsOn"]) > 0:
            pDebug("\t1> submitJob: " + job_name + " depends on " + getIDsAndNames(submitHolds)['jobnames'])
        else:
            pDebug("\t1> submitJob: " + job_name + " does not depend on other jobs" )

    # array job or single job
    if arrayJob:
        subName = job_name + "_" + str(noJobs)
        jobParams["at"] = "1"
        jobParams['lf'] = trackID + ".task"
        pDebug("\t1> submitJob: " + subName + " is an array job")
        pDebug("\t1>\tNo. tasks: " + str(noJobs))
        pDebug("\t1>\tFIRST_INDEX: " + str(taskList[0]))
        if not print_only:
            pInfo("Job: " + job_name + " - submitting an array job to " + queue + " ...")
            try:
                subOut = batchC.submit_job(
                               jobName = subName,
                               jobQueue = queue,
                               arrayProperties = { "size": noJobs },
                               jobDefinition = submitOpts["jobdef"],
                               parameters = jobParams,
                               dependsOn = submitOpts["dependsOn"],
                               containerOverrides = {
                                  "vcpus": submitOpts["vcpus"],
                                  "memory": submitOpts["memory"],
                                  "environment": submitOpts["env"]
                               },
                               retryStrategy = retryStrategy
                )
            except Exception as e:
                pError('boto3 session or client exception ' + str(e))
                sys.exit(2)
    else:
        jobParams["at"] = "0"
        jobParams['lf'] = trackID
        subName = job_name
        pDebug("\t1> submitJob: " + subName + " is a single job")
        if array_range is not None:
            pDebug("\t1> SGE_TASK_ID: " + str(taskList[0]))
        if not print_only:
            pInfo("Job: " + job_name + " - submitting single job to " + queue + " ...")
            try:
                subOut = batchC.submit_job(
                               jobName = subName,
                               jobQueue = queue,
                               jobDefinition = submitOpts["jobdef"],
                               parameters = jobParams,
                               dependsOn = submitOpts["dependsOn"],
                               containerOverrides = {
                                  "vcpus": submitOpts["vcpus"],
                                  "memory": submitOpts["memory"],
                                  "environment": submitOpts["env"]
                               },
                               retryStrategy = retryStrategy
                )
            except Exception as e:
                pError('boto3 session or client exception ' + str(e))
                sys.exit(2)
    if print_only and verbose:
        print("+++++++++  print_only verbose +++++++++++")
        print("Job: " + job_name)
        print("\tSubmit job: " + subName)
        if arrayJob:
            print("\tsubmitJob: " + subName + " is an array job")
            print("\t\tNo. tasks: " + str(noJobs))
            print("\t\tFIRST_INDEX: " + str(taskList[0]))
        elif array_range is not None:
            print("\tsubmitJob: " + subName + " is like array job but with 1 task: ")
            print("\t\tSGE_TASK_ID: " + str(taskList[0]))
        else:
            print("\tsubmitJob: " + subName + " is a single job")
        print("\tlog file: " + jobParams['lf'])
        print("\tJOB_ID: " + trackID)
        print("\tbatch queue: " + queue)
        print("\tjob definition: " + submitOpts["jobdef"])
        print("\tjob memory: " + str(submitOpts["memory"]))
        print("\tjob vcpus: " + str(submitOpts["vcpus"]))
        print("\tjob env: \n\t\t" + str(submitOpts["env"]))
        print("\tjob params: \n\t\t" + str(jobParams))
        jobid = "111-222-333-print_only-" +  subName
        subOut = {'jobName': subName, 'jobId': jobid}
        submit_id = {job_name: [jobid]}
        print("\tsubmit_id: " + str(submit_id))

    # return the "submit_id" which is a list of dictionaries
    submit_id = {job_name: [subOut['jobId']]}
    # return the job id (either from the single job or array job)
    pDebug("\t1> submitJob: " + job_name + " returning submit_id: " + str(submit_id))
    # output info file msg
    if not print_only and a_submitParams['infofile'] != None:
        with open(a_submitParams['infofile'], "a") as jifile:
            jifile.write("jobName: " + job_name + " jobQueue: " + queue + " jobId: " + subOut['jobId'] + "\n")
#   analysis log
    if a_submitParams['analysislog'] != None:
        a_submitParams['analysislog'] = lmsg

    return submit_id

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
def Summary(hdr):
    print(hdr)
    print('\tVersion: ' + fileversion)
    print('\tWorking dir: ' + str(workdir))
    print('\tCustom cluster cfg file: ' + str(cluster_file))

    print('\tAnalysis:')
    print('\t\tJob name: ' + jobname)
    print('\t\tBase cmd: ' + basecmd)
    print('\t\tPath to pipeline code: ' + apath)
    print('\t\tBase cmd args: ' + parameters)
    print('\t\tArray range: ' + str(arrayrange))
    print('\t\tNo. of cores: ' + str(nocores))
    if maxmem != None:
        print('\t\tMax memory: ' + maxmem)
    else:
        print('\t\tMax memory: specifed in cfg file')
    if queue != None:
        print('\t\tBatch queue: ' + queue)
    else:
        print('\t\tBatch queue: specified in cfg file')
    if profile != None:
        print('\t\tAWS profile: ' + profile)
    else:
        print('\t\tAWS profile: specified in cfg file')

    print('\tVerbose: ' + str(verbose))
    tbegin=time.asctime()
    print('\tTime: ' + tbegin + "\n")

def submitInteractive():

    # command line parser
    parser = ArgumentParser(description = "Helper function to submit a batch to run an analysis from analysis pipeline")
    parser.add_argument("-w", "--workdir",
                         help = "working directory (full path) [default: current working directory]")

    parser.add_argument("-j", "--jobname", help = 'Job name for tracking [default: must be specified]')
    parser.add_argument("--jobdef", help = 'Batch job definition [default: specified in config file]')
    parser.add_argument("-p", "--parameters",
                         help = 'Base cmd parameters for analysis (e.g., "-c /usr/local/analysis_pipeline/R/ld_pruning.R cfgfile.cfg --version xxx")')
    parser.add_argument("--apath", help = "analysis pipeline path [default: " + defApath + "]")
    parser.add_argument("-b", "--basecmd", help = "Base command for analysis [default: " + defBaseCmd + "]")

    parser.add_argument("--cluster_file", help = "custom batch config file [default: None]")
    parser.add_argument("--arrayrange", help = "job array range (e.g., 1-22) [default: None]")
    parser.add_argument("--nocores", help = "Number of cores [default: single core]")
    parser.add_argument("-M", "--maxmem", help = "Maximum memory (MB) [default: specified in config file]")
    parser.add_argument("-q", "--queue",  help = "batch queue [default: specified in config file]")
    parser.add_argument("-P", "--profile", help = "AWS profile [default: specified in batch config file]")

    parser.add_argument("-V", "--verbose", action="store_true", default = False,
                         help = "Turn on verbose output [default: False]")
    parser.add_argument("-S", "--summary", action="store_true", default = False,
                         help = "Print summary prior to executing [default: False]")
    parser.add_argument("--print_only", action="store_true", default = False,
                         help = "Test without executing [default: False]")
    parser.add_argument("--version", action="store_true", default = False,
                         help = "Print version of " + __file__)

    global jobname
    global jobdef
    global basecmd
    global parameters
    global workdir
    global apath
    global cluster_file
    global nocores
    global maxmem
    global queue
    global profile
    global verbose
    global summary
    global print_only
    global arrayrange
    global version
    args = parser.parse_args()
    jobname = args.jobname
    jobdef = args.jobdef
    basecmd = args.basecmd
    parameters = args.parameters
    workdir = args.workdir
    apath = args.apath
    cluster_file = args.cluster_file
    nocores = args.nocores
    maxmem = args.maxmem
    queue = args.queue
    profile = args.profile
    verbose = args.verbose
    summary = args.summary
    print_only = args.print_only
    arrayrange = args.arrayrange
    version = args.version

    if args.version:
        print(__file__ + " version: " + fileversion)
        sys.exit()
    # job name
    if jobname == None:
        pError("jobname option (-j/--jobname) has not be specified")
        sys.exit(2)
    # check on analysis and parameters (required params)
    if parameters == None:
        pError("Base cmd's parameters for analysis (-p argument) must be specified")
        sys.exit(2)

    if cluster_file != None:
        if not os.path.isfile(cluster_file):
            pError("Cluster config file " + cluster_file + " does not exist")
            sys.exit(2)

    # get the cluster configuration
    cfgversion = "3.2"
    if apath == None:
        cpath = defApath
    else:
        cpath = apath
    stdcfgfile = cpath + "/aws_batch_cfg.json"
    clustercfg = getClusterCfg(stdcfgfile, cluster_file, cfgversion)
    if summary:
        Summary("Summary of " + __file__)

    # submit (with print_only flag to not really submit)
    global subParams
    subParams['profile'] = profile
    subParams['cluster_file']= cluster_file
    subParams['clustercfg'] = clustercfg
    subParams['jobname'] = jobname
    subParams['cmd'] = basecmd
    subParams['args'] = parameters
    subParams['holdid'] = None
    subParams['array_range'] = arrayrange
    subParams['request_cores'] = nocores
    subParams['print_only'] = print_only
    subParams['verbose'] = verbose
    subParams['maxmem'] = maxmem
    subParams['workdir'] = workdir
    subParams['queue'] = queue
    subParams['jobdef'] = jobdef
    subParams['apath'] = apath

    jobid = submitjob(subParams)
    if not print_only:
        pInfo('Job ' + str(jobid) + ' submitted.')

if __name__ == "__main__":
    submitInteractive()
