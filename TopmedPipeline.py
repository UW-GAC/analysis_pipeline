"""Utility functions for TOPMed pipeline"""

import sys
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
