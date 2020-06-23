#! /usr/bin/env python3
import     time
import     csv
import     sys
import  os.path
import  os
import  subprocess
from    argparse import ArgumentParser
from     datetime import datetime, timedelta

# init globals
version='1.4'
msgErrPrefix='>>> Error: '
msgInfoPrefix='>>> Info: '
debugPrefix='>>> Debug: '

# def functions
def pInfo(msg):
    print(msgInfoPrefix+msg)

def pError(msg):
    print(msgErrPrefix+msg)

def pDebug(msg):
    if debug:
        print(debugPrefix+msg)
def Summary(hdr):
    pInfo(hdr)
    pInfo( '\tVersion: ' + version)
    pInfo( '\tR driver file: ' + rdriver )
    pInfo( '\tR driver args:' + rargs)
    if logfile != "":
        pInfo( '\tLog file:' + fullLog)
    pInfo( '\tSystem parameters:' )
    pInfo( '\t\tWorking directory: ' + workdir )
    pInfo( '\t\tData root folder: '+ dataroot )
    if nomount:
        pInfo("\t\tNot mounting network storage (should be mounted via docker run -v command)")
    else:
        pInfo( '\t\tMount command: '+ mount )
        pInfo( '\t\tMount timeout: '+ tmo)
    if debug:
        pInfo( '\tDebug: True' )
    else:
        pInfo( '\tDebug: False' )
    if arrayType:
        pInfo( '\tArray type: True' )
        pInfo( '\tArray index: ' + arrayIndex )
        pInfo( '\tFirst index: ' + firstSegIndex)
        pInfo( '\tSGE Task ID: ' +  os.environ['SGE_TASK_ID'])
    else:
        pInfo( '\tArray type: False')
    tbegin=time.asctime()
    pInfo( '\tTime: ' + tbegin + "\n" )


# default names
defaultRdriver = "/usr/local/analysis_pipeline/runRscript.sh"
defaultMount = "mount -t nfs4 -o vers=4.1 172.255.44.97:/ /projects"
defaultDataroot = "/projects"
defaultWorkdir = defaultDataroot + "/batch"
defaultRargs = ""
defaultDebug = 0
defaultPO = 0
defaultArrayType = 0

# command line parser
parser = ArgumentParser( description = "docker script to run tm analysis pipeline R code via R control file (e.g., runRscript.sh)" )
parser.add_argument( "-w", "--workdir", default = defaultWorkdir,
                     help = "full path of working directory (where R runs) [default: " + defaultWorkdir + "]" )
parser.add_argument( "-d", "--dataroot", default = defaultDataroot,
                     help = "external data root folder [default: " + defaultDataroot + "]" )
parser.add_argument( "-n", "--nomount", help = "do not mount nfs/cifs disk [default: false]",
                     action = "store_true")
parser.add_argument("--rdriver", default = defaultRdriver,
                     help = "full path of pipeline R driver bash file [default: " + defaultRdriver + "]" )
parser.add_argument("-m", "--mount", default = defaultMount,
                     help = "nfs mount of data volume [default: "+defaultMount+"]" )
parser.add_argument("--rargs", default = defaultRargs,
                     help = "R driver arguments" )
parser.add_argument( "-D", "--Debug", type = int, default = defaultDebug,
                     help = "Turn on debug output [default: " + str(defaultDebug) + "]" )
parser.add_argument( "-p", "--printonly", type = int, default = defaultPO,
                     help = "Print summary without executing [default: " + str(defaultPO) + "]" )
parser.add_argument( "-l", "--logfile", help = "log filename", default = "" )
parser.add_argument( "-t", "--timeout", help = "mount timeout", default = "10.0" )
parser.add_argument( "-a", "--arraytype", type = int, default = defaultArrayType,
                     help = "Batch job is array type [default: " + str(defaultArrayType) + " ]" )
parser.add_argument("--version", action="store_true", default = False,
                    help = "Print version of " + __file__ )


args = parser.parse_args()
# set result of arg parse_args
rdriver = args.rdriver
workdir = args.workdir
dataroot = args.dataroot
rargs = args.rargs
mount = args.mount
logfile = args.logfile
tmo = args.timeout

debug = args.Debug
po = args.printonly
nomount = args.nomount
arrayType = args.arraytype
# version
if args.version:
    print(__file__ + " version: " + version)
    sys.exit()

# handle array type
if arrayType:
    # get the batch array index, the first segment index and set SGE_TASK_ID
    # which corresponds to a segment index in a chromosome
    echeck = 'AWS_BATCH_JOB_ARRAY_INDEX'
    if echeck in os.environ:
        arrayIndex = os.environ[echeck]
    else:
        pError("Required array job env " + echeck + " not found")
        sys.exit(2)
    echeck = 'FIRST_INDEX'
    if echeck in os.environ:
        firstSegIndex = os.environ[echeck]
    else:
        pError("Required first index env " + echeck + " not found")
        sys.exit(2)
    taskID = str(int(arrayIndex) + int(firstSegIndex))
    os.environ['SGE_TASK_ID'] = taskID

# check for logile; if so, make it a full path to working directory
logExt = ".log"
if logfile != "":
    if arrayType:
        logfile = logfile + "_" + taskID + logExt
    else:
        logfile = logfile + logExt
    fullLog = workdir + "/" + logfile

# summarize and check for required params
Summary("Summary of " + __file__)

if po:
    pInfo( "Exiting without executing." )
    sys.exit()

# check if the mount point (last arg in mount command) exists; if not create it
if not os.path.isdir( mount.split()[-1] ):
    os.mkdir(mount.split()[-1])
# mount
if not nomount:
    pDebug( "mount tmo: " + tmo + " mount command: " + mount )
    sys.stdout.flush()
    mtmo = "timeout " + tmo + " " + mount
    process = subprocess.Popen(mtmo, shell=True, stdout=subprocess.PIPE)
    status = process.wait()
    pipe = process.stdout
    msg = pipe.readline()
    if status == 32:
        pInfo("Warning: mount volume already mounted.")
    elif status:
        if status == 124:
            msg = "mount timed out"
        pError("Mount error: " + msg )
        sys.exit(2)
else:
    pDebug( "Not mounting network storage" )

# check for files and folders exists
if not os.path.isfile( rdriver ):
    pError( "R control file " + rdriver + " does not exist" )
    sys.exit(2)
if not os.path.isdir( workdir ):
    pError( "Work directory " + workdir + " does not exist" )
    sys.exit(2)
if not os.path.isdir( dataroot ):
    pError( "Data root directory " + dataroot + " does not exist" )
    sys.exit(2)
# change working directory
os.chdir(workdir)
pInfo( "Changed working directory to " + workdir )
# execute the R control file
cmd = rdriver + ' ' + rargs
pDebug( "Executing " + cmd )
sys.stdout.flush()
# redirect stdout to logile
if logfile != "":
    sout = sys.stdout
    serr = sys.stderr
    flog = open ( fullLog, 'w' )
    sys.stderr = sys.stdout = flog
process = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr, shell=True)
status = process.wait()
# redirect stdout back
if logfile != "":
    sys.stdout = sout
    sys.stderr = serr
if status:
    pError( "Executing R driver file failed (" + str(status) + ")" )
    sys.exit(2)
else:
    pInfo( "Executing R driver file completed without errors.")
