#!/usr/bin/python
import     time
import     csv
import     sys
import     os.path
import     os
import     subprocess
from       argparse import ArgumentParser
from       datetime import datetime, timedelta
from       shutil import copyfile


# init globals
version='1.0'
msgErrPrefix='>>> Error: '
msgInfoPrefix='>>> Info: '
debugPrefix='>>> Debug: '

# def functions
def pInfo(msg):
    print msgInfoPrefix+msg

def pError(msg):
    print msgErrPrefix+msg

def pDebug(msg):
    if debug:
        print debugPrefix+msg
def Summary(hdr):
    pInfo(hdr)
    pInfo( '\tVersion: ' + version)
    pInfo( '\tWork directory: ' + workDir )
    pInfo( '\tPipeline command:' + pipelineCmd)
    pInfo( '\tPipeline args:' + pipelineArgs)
    pInfo( '\tAWS cli security file:' + awsSecurity)
    if logfile != "":
        pInfo( '\tLog file:' + fullLog)
    if debug:
        pInfo( '\tDebug: True' )
    else:
        pInfo( '\tDebug: False' )
    tbegin=time.asctime()
    pInfo( '\tTime: ' + tbegin + "\n" )


# default names
defaultDebug = 0
defaultPO = 0

# command line parser
parser = ArgumentParser( description = "docker script to run tm analysis pipeline R code via R control file (e.g., runRscript.sh)" )
parser.add_argument( "-w", "--workdir",
                     help = "full path of working directory (where pipeline jobs are submitted)" )
parser.add_argument("-c","--command",
                     help = "pipeline command" )
parser.add_argument("-a","--args",
                     help = "pipeline command arguments" )
parser.add_argument("-s", "--securityfile",
                     help = "full path of aws security credentials file" )
parser.add_argument( "-D", "--debug", type = int, default = defaultDebug,
                     help = "Turn on debug output [default: " + str(defaultDebug) + "]" )
parser.add_argument( "-P", "--printonly", type = int, default = defaultPO,
                     help = "Print summary without executing [default: " + str(defaultPO) + "]" )
parser.add_argument( "-l", "--logfile", help = "log filename" )
parser.add_argument("--version", action="store_true", default = False,
                    help = "Print version of " + __file__ )


args = parser.parse_args()
# set result of arg parse_args
workDir = args.workdir
pipelineCmd = args.command
pipelineArgs = args.args
awsSecurity = args.securityfile
logfile = args.logfile
debug = args.Debug
po = args.printonly
# version
if args.version:
    print(__file__ + " version: " + version)
    sys.exit()

# check for logile; if so, make it a full path to working directory
logExt = ".log"
if logfile != "":
    if arrayType:
        logfile = logfile + "_" + taskID + logExt
    else:
        logfile = logfile + logExt
    fullLog = workDir + "/" + logfile

# summarize and check for required params
Summary("Summary of " + __file__)

if po:
    pInfo( "Exiting without executing." )
    sys.exit()

# check for files and folders exists
if not os.path.isfile( rdriver ):
    pError( "R control file " + rdriver + " does not exist" )
    sys.exit(2)
if not os.path.isdir( workDir ):
    pError( "Work directory " + workDir + " does not exist" )
    sys.exit(2)

# change working directory
os.chdir(workDir)
pInfo( "CD to " + workDir )

# copy to security file
pInfo( "Coping " + awsSecurity " to ~/.aws " )
cmd = 'copy ' + awsSecurity ' ~/.aws'
copyfile(awsSecurity, '~/.aws')
# execute the pipeline command
cmd = pipelineCnmd + ' ' + pipelineArgs
pInfo( "Executing " + cmd )
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
    pError( "Executing pipeline command failed (" + str(status) + ")" )
    sys.exit(2)
else:
    pInfo( "Executing pipeline command completed without errors.")
