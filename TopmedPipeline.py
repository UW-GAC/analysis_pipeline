"""Utility functions for TOPMed pipeline"""

import sys
import csv
import subprocess
from copy import deepcopy
import getpass
import time

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


def getOptions(file):
    """Read a file in .sge_request format and return a dictionary"""
    opts = {}
    if file is not None:
        f = open(file, 'r')
        for line in f:
            vals = line.split()
            if len(vals) == 1:
                vals.append('')
            (key, value) = vals
            opts[key] = value
        f.close()
    return opts



# parent class to represent a compute cluster environment
class Cluster(object):
    """ """
    # constructor
    def __init__(self, submit_cmd, options=None):
        if options is None:
            options = {}

        self.submit_cmd = submit_cmd
        self.options = options


    def submitJob(self, cmd, args=None, opts=None, verbose=True, print_only=False):
        """ Submit a job to the cluster and return job id"""

        if args is None:
            args = []
        argStr = " ".join(args)

        # override any stored options with argument
        if opts is None:
            opts = {}
        options = deepcopy(self.options)
        options.update(opts)
        optStr = dictToString(options)

        sub_cmd = " ".join([self.submit_cmd, optStr, cmd, argStr])

        if print_only:
            print sub_cmd
            return "000000"

        process = subprocess.Popen(sub_cmd, shell=True, stdout=subprocess.PIPE)
        pipe = process.stdout
        sub_out = pipe.readline()
        jobid = sub_out.split()[2]

        if verbose:
            print sub_out

        return jobid

class AWS_Batch(Cluster):

    def __init__(self, options=None):
        user = getpass.getuser()
        wdkey = "wd"
        if options is not None:
            if wdkey in options:
                defwd = options[wdkey]
            else:
                defwd = os.getcwd()
        self.defParams = { wdkey: defwd,
                           "dr": "/projects",
                           "rd": "/usr/local/analysis_pipeline/runRscript.sh",
                           "ra": "",
                           "db": "1",
                           "po": "0",
                           "lf": "",
                           "mt": "mount -t nfs4 -o vers=4.1 172.31.58.33:/ /projects"
                         }
        self.defQueue = "LowPriority"
        self.defJobdef = "topmed_general"
        self.defVcpus =  2
        self.defMemory =  8000
        self.defRegion = "us-east-1"
        self.defDependsOn = []
        self.defEnv = []
        self.memkey = "-m"
# create the batch client
        self.batchC = boto3.client('batch',region_name=self.defRegion)

    def memoryOptions(self, job_name):
        vmem = {"find_unrelated":3.5,
                "ld_pruning":11000,
                "combine_variants":1000,
                "pca_byrel":4000,
                "pca_plots":1000,
                "pca_corr":7000,
                "pca_corr_plots":32000,
                "pcrelate":8000,
                "null_model":36000,
                "aggregate_list":6000,
                "assoc_single":24000,
                "assoc_aggregate":24000,
                "assoc_window":24000,
                "assoc_combine":6000,
                "assoc_plots":4000}
        a =[ akey for akey in vmem.keys() if job_name.startswith( akey ) ]
        if len(a):
            opts = { self.memkey: vmem[a[0]] }
        else:
            opts = {}
        return opts

    def submitJob(self, job_name, cmd, args=None, holdid=None, array_range=None, opts=None, print_only=False, **kwargs):

# job def/name
        self.defJobName = job_name

# process the args/kwargs to go the driver (rd); the driver arguments (e.g., -s rcode cfg --chr cn)
        self.defParams['rd'] = cmd
        if args is None:
            args = []
        self.defParams['ra'] = " ".join(args)

# assign the log file
        t = str(int(time.time()))
        lfile = job_name + "_" + t + ".log"
        self.defParams['lf'] = lfile

# print_only
        if print_only:
            self.defParams['po'] = "1"

# get memory option
        if opts is None:
            opts = {}
        else:
            if self.memkey in opts:
                self.defMemory = opts[self.memkey]
# holdid is a list of a list of jid dictionaries
        if holdid is not None and holdid[0] != []:
            self.defDependsOn = holdid[0]

# process kwargs to get

# if we're doing a job array equivalent submit multiple jobs; else just submit one job
        idsonly = []
        jids = []
        if not print_only:
            if array_range is not None:
                lf = self.defParams['lf']
                air = [ int(i) for i in array_range.split( '-' ) ]
                for t in range( air[0], air[1]+1 ):
                    # add an environmnent for equiv to taskid in sge
                    self.defEnv = [ { "name": "SGE_TASK_ID",
                                      "value": str(t) } ]
                    self.defParams['lf'] = lf + '.' + str(t)

                    subOut = self.batchC.submit_job(
                       jobName = self.defJobName,
                       jobQueue = self.defQueue,
                       jobDefinition = self.defJobdef,
                       parameters = self.defParams,
                       dependsOn = self.defDependsOn,
                       containerOverrides = {
                          "vcpus": self.defVcpus,
                          "memory": self.defMemory,
                          "environment": self.defEnv
                       }
                    )
                    # append the jid dictionary to list
                    jids.append( { "jobId": subOut["jobId"] } )
                    idsonly.append( subOut["jobId"] )
            else:
                subOut = self.batchC.submit_job(
                   jobName = self.defJobName,
                   jobQueue = self.defQueue,
                   jobDefinition = self.defJobdef,
                   parameters = self.defParams,
                   dependsOn = self.defDependsOn,
                   containerOverrides = {
                      "vcpus": self.defVcpus,
                      "memory": self.defMemory,
                      "environment": self.defEnv
                   }
                )
                # append the jid dictionary to list
                jids.append( { "jobId": subOut["jobId"] } )
                idsonly.append( subOut["jobId"] )
            print( "Job Name: " + self.defJobName + " Job IDs: " + str(idsonly) )
        else:
            batchCmd = "\n\tjobName = " + self.defJobName
            batchCmd = batchCmd + ", jobQueue = " + self.defQueue
            batchCmd = batchCmd + ", jobDef = " + self.defJobdef
            batchCmd = batchCmd + ", memory = " + str(self.defMemory)
            batchCmd = batchCmd + ", vcpus = " + str(self.defVcpus)
            batchCmd = batchCmd + "\n\tparams = " + str(self.defParams)
            if array_range is not None:
                ct = "array job " + array_range
            else:
                ct = "single job"
            if holdid is None or holdid[0] == []:
                ht = "no hold"
            else:
                ht = holdid[0]

            print( "\njob name: " + self.defJobName + " job type: " + ct + "\tht: " + str(ht) + "\nbatch.submit_job : " + batchCmd)
            jids = [ { 'jobId': 'print_only-1/' + self.defJobName },
                     { 'jobId': 'print_only-2/' + self.defJobName } ]

# return a list of jid dictionaries [ {'jobId': string_of_jobid}, ...]
        return jids


class SGE_Cluster(Cluster):

    def __init__(self, options=None, parallel_environment="local"):
        super(SGE_Cluster, self).__init__(submit_cmd="qsub", options=options)
        self.pe = parallel_environment


    def submitJob(self, job_name, holdid=None, array_range=None, request_cores=None, email=None, opts=None, **kwargs):

        if opts is None:
            opts = {}

        opts["-N"] = job_name

        if holdid is not None and holdid != []:
            if isinstance(holdid, str):
                holdid = [holdid]
            opts["-hold_jid"] =  ",".join(holdid)

        if array_range is not None:
            opts["-t"] = array_range

        if request_cores is not None:
            opts["-pe"] = self.pe + " " + request_cores

        if email is not None:
            opts["-m"] = "e"
            opts["-M"] = email

        jobid = super(SGE_Cluster, self).submitJob(opts=opts, **kwargs)

        if array_range is not None:
            jobid = jobid.split(".")[0]

        return jobid


    def memoryOptions(self, job_name):
        # requesting memory causes problems on the UW cluster
        opts = {}
        return opts



class UW_Cluster(SGE_Cluster):

    def __init__(self, options=None):
        defaults = {"-cwd":"",
                    "-j":"y",
                    "-q":"olga.q",
                    "-S":"/bin/bash",
                    "-v":"R_LIBS=/projects/resources/gactools/R_packages/library,PATH=/projects/resources/software/apps/bin:$PATH"}

        if options is None:
            options = {}
        defaults.update(options)

        super(UW_Cluster, self).__init__(options=defaults, parallel_environment="local")



class AWS_Cluster(SGE_Cluster):

    def __init__(self, options=None):
        defaults = {"-cwd":"",
                    "-j":"y",
                    "-q":"all.q",
                    "-S":"/bin/bash",
                    "-v":"R_LIBS=/projects/resources/gactools/R_packages/library,PATH=/projects/resources/software/apps/bin:$PATH"}

        if options is None:
            options = {}
        defaults.update(options)

        super(AWS_Cluster, self).__init__(options=defaults, parallel_environment="local")


    def submitJob(self, *args, **kwargs):
        # currently, no email on aws
        kwargs["email"] = None
        jobid = super(AWS_Cluster, self).submitJob(*args, **kwargs)
        return jobid


    def memoryOptions(self, job_name):
        vmem = {"find_unrelated":3.5,
                "ld_pruning":11,
                "combine_variants":1,
                "pca_byrel":4,
                "pca_plots":1,
                "pca_corr":5.5,
                "pca_corr_plots":132,
                "null_model":6,
                "aggregate_list":6,
                "assoc":4.5,
                "assoc_combine":3,
                "assoc_plots":1.5}

        option = "-l"
        resource = "h_vmem={0}G".format(vmem[job_name])
        opts = {}
        opts[option] = resource
        return opts



class ClusterFactory(object):

    @staticmethod
    def createCluster(cluster_type, *args, **kwargs):
        if cluster_type == "uw":
            return UW_Cluster(*args, **kwargs)
        elif cluster_type == "aws":
            return AWS_Cluster(*args, **kwargs)
        else:
            raise Exception("unknown cluster type: " + cluster_type + "!")
