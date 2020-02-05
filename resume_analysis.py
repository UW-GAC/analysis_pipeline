#! /usr/bin/env python3
# Support resume on slurm.
# Monitors the jobs/tasks of analysis.  If enabled, resume will do the following:
#   1.  Determines if this job is an array job (cluster dependent)
#   2.  Creates a file for each job that completes and runs without errors
#   3.  If a completed/success file for a job exists, then the this script just returns
#   4.  Else this scrip runs the job (typically an R job) and if it completes without
#       it then touches the complete file
#
# The name of completed file associated with the job is:
#       <jobname>_<task_id>.completed (for an array job)
#       <jobname>.completed (for a single job)
# The script will place the job id into the completed file.
#
# There are 4 or more arguments:
#       <resumescript> <cluster_type> <jobname> <job run cmd>  <job run cmd arguments>
# This script basically parses out the jobrun command and all it's arguments, and then
# calls popen to execute syncrhonously.  Once completed, this script checks the status and
# touches the completed file.
#
# and the sc
from __future__ import print_function
import      sys
import      errno
import      os
import      time
import      subprocess
import      glob

import      port_popen

# init globals
fileversion = '1.0'
msgErrPrefix = '>>> Error: '
msgInfoPrefix = '>>> Info: '
debugPrefix = '>>> Vebose: '
verbose = False

class JobScheduler(object):
    def __init__(self):
        print("JobScheduler init ...")
        self.jobInfo = {"arrayjob":None, "jobid": None, "taskid": None}
        self.resumeDir = "./resume/"
    def getJobInfo(self):
        return self.jobInfo
    def updateJobInfo(self, jinfo):
        self.jobInfo = jinfo
    def setEnvInfo(self, envInfo):
        for key in list(envInfo.keys()):
            envInfo[key] = os.getenv(key)


class SGE(JobScheduler):
    def __init__(self):
        super(SGE,self).__init__()
        self.envInfo = {
                "JOB_ID": None,
                "SGE_TASK_ID": None
        }
        self.setJobInfo()

    def setJobInfo(self):
        jobInfo = super(SGE, self).getJobInfo()
        super(SGE, self).setEnvInfo(self.envInfo)
        # if array job, use the array job id; else job id
        # check if array_job
        if self.envInfo["SGE_TASK_ID"] != "undefined":
            jobInfo["arrayjob"] = True
            jobInfo["jobid"] = self.envInfo["JOB_ID"]
            jobInfo["taskid"] = self.envInfo["SGE_TASK_ID"]
        else:
            jobInfo["arrayjob"] = False
            if self.envInfo["JOB_ID"] != None:
                jobInfo["jobid"] =  self.envInfo["JOB_ID"]
            else:
                jobInfo["jobid"] = "NOJOBID"
        super(SGE, self).updateJobInfo(jobInfo)

    def createFilename(self, jName, ctag=None):
        jInfo = super(SGE, self).getJobInfo()
        if ctag == None:
            pfile = jName + ".o" + jInfo["jobid"]
            if jInfo["arrayjob"]:
                pfile += "." + jInfo["taskid"]
        else:
            if jInfo["arrayjob"]:
                pfile = ctag + "_" + jName + "_" + jInfo["taskid"]
            else:
                pfile = ctag + "_" + jName

        return pfile


class BATCH(JobScheduler):
    def __init__(self):
        super(BATCH,self).__init__()
        self.envInfo = {
            "AWS_BATCH_JOB_ARRAY_INDEX": None,
            "AWS_BATCH_JOB_ID": None
        }
        self.setJobInfo()

    def setJobInfo(self):
        jobInfo = super(BATCH, self).getJobInfo()
        super(BATCH, self).setEnvInfo(self.envInfo)
        # if array job, use the array job id; else job id
        # check if array_job
        if self.envInfo["AWS_BATCH_JOB_ARRAY_INDEX"] != None:
            jobInfo["arrayjob"] = True
            jobInfo["jobid"] = self.envInfo["AWS_BATCH_JOB_ID"]
            jobInfo["taskid"] = self.envInfo["AWS_BATCH_JOB_ARRAY_INDEX"]
        else:
            jobInfo["arrayjob"] = False
            if self.envInfo["AWS_BATCH_JOB_ID"] != None:
                jobInfo["jobid"] =  self.envInfo["AWS_BATCH_JOB_ID"]
            else:
                jobInfo["jobid"] = "NOJOBID"
        super(BATCH, self).updateJobInfo(jobInfo)

    def createFilename(self, jName, ctag=None):
        jInfo = super(BATCH, self).getJobInfo()
        if ctag == None:
            pfile = jName + "_" + jInfo["jobid"] + ".log"
            if jInfo["arrayjob"]:
                pfile += "_" + jInfo["taskid"] + ".log"
        else:
            if jInfo["arrayjob"]:
                pfile = ctag + "_" + jName + "_" + jInfo["taskid"]
            else:
                pfile = ctag + "_" + jName
        return pfile

class SLURM(JobScheduler):
    def __init__(self):
        super(SLURM,self).__init__()
        self.envInfo = {
            "SLURM_ARRAY_TASK_ID": None,
            "SLURM_ARRAY_JOB_ID": None,
            "SLURM_JOB_ID": None
        }
        self.setJobInfo()

    def setJobInfo(self):
        jobInfo = super(SLURM, self).getJobInfo()
        super(SLURM, self).setEnvInfo(self.envInfo)
        # if array job, use the array job id; else job id
        # check if array_job
        if self.envInfo["SLURM_ARRAY_JOB_ID"] != None:
            jobInfo["arrayjob"] = True
            jobInfo["jobid"] = self.envInfo["SLURM_ARRAY_JOB_ID"]
            jobInfo["taskid"] = self.envInfo["SLURM_ARRAY_TASK_ID"]
        else:
            jobInfo["arrayjob"] = False
            if self.envInfo["SLURM_JOB_ID"] != None:
                jobInfo["jobid"] =  self.envInfo["SLURM_JOB_ID"]
            else:
                jobInfo["jobid"] = "NOJOBID"
        super(SLURM, self).updateJobInfo(jobInfo)

    def createFilename(self, jName, ctag=None):
        jInfo = super(SLURM, self).getJobInfo()
        if ctag == None:
            pfile = jName + "_" + jInfo["jobid"] + ".log"
            if jInfo["arrayjob"]:
                pfile += "_" + jInfo["taskid"] + ".log"
        else:
            if jInfo["arrayjob"]:
                pfile = ctag + "_" + jName + "_" + jInfo["taskid"]
            else:
                pfile = ctag + "_" + jName
        return pfile

def pInfo_file(msg, fhandle):
    tmsg=time.asctime()
    print(msgInfoPrefix+tmsg+": "+msg, file=fhandle)

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
    print("\tVersion: " + fileversion)

# get the args
minArgs = 4
if (len(sys.argv) < minArgs):
    pError("There must be at least " + str(minArgs) + " arguments provided to this script.")
    sys.exit(2)
sName = sys.argv[0]
cName = sys.argv[1]
jName = sys.argv[2]
# set the cmd
jCmd = sys.argv[3]
if (len(sys.argv) > minArgs):
    # must iterate over the arguments in case an argument contains a space
    jArgs = ""
    for a in sys.argv[minArgs:len(sys.argv)]:
        if a.find(' ') == -1:
            jArgs += ' ' + a
        else:
            jArgs += " '" + a + "'"
    jCmd += ' ' + jArgs
# check and instantiate cluster class
clusters = ['SLURM', 'BATCH', 'SGE']
if cName.upper() not in clusters:
    pError(sName + "- Error: invalid cluster type " + cName)
    sys.exit(2)
clusterclass = globals()[cName.upper()]
thecluster = clusterclass()
pInfo("Script: " + sName)
pInfo("Job: " + jName)
pInfo("Cmd: " + jCmd)
# get the job info
jInfo = thecluster.getJobInfo()
jId = jInfo["jobid"]
# set the file name
ctag = "completed"
cfname = thecluster.resumeDir + thecluster.createFilename(jName, ctag)
lfname = thecluster.createFilename(jName)

status = 0

if os.path.isfile(cfname):
    pInfo("Job " + jName + " already completed (" + cfname +")")
else:
    pInfo("Calling Popen to execute: \n\t" + jCmd)
    (pmsg, status) = port_popen.popen_stdout(jCmd, logfile=lfname, jobname=jName)
    # if status == 0, write a complete file unless it's the post_analysis
    # open the file and write out job info (name, id, task)
    if status == 0 and "post_analysis" not in jName:
        with open(cfname, mode='w') as fhandle:
            msgC =  "job: " + jName + ", job id: " + jId
            if jInfo["arrayjob"]:
                pInfo_file(msgC + " task: " + jInfo["taskid"] + " completed", fhandle)
            else:
                pInfo_file(msgC + " completed", fhandle)
    else:
        pError("Status: " + str(status) + " See " + lfname + " for details")

sys.exit(status)
