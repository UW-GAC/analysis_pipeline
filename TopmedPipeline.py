"""Utility functions for TOPMed pipeline"""

import csv
import subprocess



def readConfig(file):
    """Read a pipeline config file.

    Usage: 
    config = readConfig(file)

    Arguments: 
    file - name of config file to read

    Returns: 
    dictionary with config values
    """

    config = dict()
    f = open(file, 'r')
    reader = csv.reader(f, delimiter=' ', quotechar='"', skipinitialspace=True)
    for line in reader:
        if line[0][0] == "#":
            continue

        if len(line) == 2:
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



def submitJob(job, cmd, args, queue, holdid=None, arrayRange=None, requestCores=None,
              email=None, qsubOptions="", verbose=True, printOnly=False):
    """Sumbit a pipeline job.

    Usage: 
    jobid = submitJob(job, cmd, args, queue, holdid=None, arrayRange=None, requestCores=None, email=None, qsubOptions="", verbose=True, printOnly=False)

    Arguments:
    job - name of job
    cmd - command to execute
    args - list of arguments to cmd
    queue - compute cluster queue to submit job into
    holdid - list of job ids that must be complete before this job is run
    arrayRange - specified, range of array jobs to pass to qsub (ie 1-23)
    requestCores - number of cores to request; either a number (e.g, 1) or a range of numbers (e.g., 1-4)
    email - email address to notify when job is complete
    qsubOptions - additional options to pass to qsub
    verbose - Print out stdout from qsub?
    printOnly - print qsub command without submitting

    Returns:
    id of the submitted job

    """

    nameStr = "-N " + job
    
    queueStr = "-q " + queue
    
    if holdid is not None and holdid != []:
        if isinstance(holdid, str):
            holdid = [holdid]
        holdStr = "-hold_jid " + ",".join(holdid)
    else:
        holdStr = ""

    if requestCores is not None:
        coreStr = "-pe local " + requestCores
    else:
        coreStr = ""
        
    if arrayRange is not None:
        arrayStr = "-t " + arrayRange
    else:
        arrayStr = ""
        
    if email is not None:
        emailStr = "-m e -M " + email
    else:
        emailStr = ""

    argStr = " ".join(args)

    qsub = "qsub -S /bin/bash %s %s %s %s %s %s %s %s %s" % (qsubOptions, nameStr, arrayStr, holdStr, coreStr, queueStr, emailStr, cmd, argStr)
    
    if printOnly:
        print qsub
        return "000000"
    
    process = subprocess.Popen(qsub, shell=True, stdout=subprocess.PIPE)
    pipe = process.stdout
    qsubout = pipe.readline()
    jobid = qsubout.split()[2]
    
    if verbose:
        print qsubout
    
    return jobid
