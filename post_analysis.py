#! /usr/bin/env python
import     time
from       argparse import ArgumentParser
import     datetime
import     fnmatch
import     os
import     sys
import     glob

# check for number
def isnumber(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# parse input
parser = ArgumentParser( description = "script for post analysis stuff" )
parser.add_argument( "-a", "--analysis",
                     help = "Name of analysis" )
parser.add_argument( "-l", "--logfile",
                     help = "Name of analysis log file" )
parser.add_argument( "-s", "--starttime",
                     help = "start time (in datetime format) of analysis ")

args = parser.parse_args()
analysis = args.analysis
starttime = args.starttime
logfile = args.logfile
# strip of any leading/trailing quotes
starttime = starttime.replace("_", " ")
# check for any errors in the log and resume files
foundErr = False
fileError = []
rfiles = glob.glob("resume*")
if len(rfiles):
    for file in rfiles:
        with open(file) as f:
            if ">>> Error" in f.read():
                fileError.append(file)
                foundErr = True;
                break;
ofiles = glob.glob("*.o*")
if len(ofiles):
    for file in ofiles:
        with open(file) as f:
            if "Error" in f.read():
                fileError.append(file)
                foundErr = True;
                break;
lfiles = glob.glob("*.log")
if len(lfiles):
    for file in lfiles:
        with open(file) as f:
            if "Error" in f.read():
                fileError.append(file)
                foundErr = True;
                break;
ffiles = glob.glob("fail*")
if len(ffiles):
    fileError.append(",".join(ffiles))
    foundErr = True;
if foundErr:
    print("post_analysis not executed because analysis error(s) has been found in " + str(fileError))
    sys.exit(2)
# end time (utc)
tFmt = "%a, %d %b %Y %I:%M:%S %p"
endtime = datetime.datetime.utcnow().strftime(tFmt)
deltaTime = (datetime.datetime.strptime(endtime, tFmt) -
             datetime.datetime.strptime(starttime, tFmt)).total_seconds()
td_hrs = deltaTime/60./60.
# append to log file
with open(logfile, "a") as afile:
    afile.write(analysis + " end time: " + endtime + "\n")
    afile.write(analysis + " elapse time(hrs): " + str(td_hrs) + "\n")
# get cost (for batch and slurm)
totCost = 0.0;
for file in os.listdir('.'):
    if file.endswith(".log"):
        with open(file) as thefile:
            for line in thefile:
                if line.startswith(">>> Info") and line.lower().find("cost"):
                    ll = line.split("$")
                    if len(ll) == 2:
                        if isnumber(ll[1]):
                            totCost += float(ll[1])
if totCost != 0.0:
    with open(logfile, "a") as afile:
        afile.write(analysis + " total estimated cost: $" + str(round(totCost,2)) + "\n")

# cleanup
errcnt = 0
# delete all the resume* and complete*
rfiles = glob.glob("resume*")
if len(rfiles):
    for file in rfiles:
        os.remove(file)

cfiles = glob.glob("completed*")
if len(cfiles):
    for file in cfiles:
        os.remove(file)

for file in os.listdir('.'):
    try:
        if fnmatch.fnmatch(file, '*.log') or fnmatch.fnmatch(file,'*.trace') or \
           fnmatch.fnmatch(file, '*.o*') or fnmatch.fnmatch(file,'*.po*'):
                os.rename('./'+file,'./log/'+file)
        elif fnmatch.fnmatch(file, '*report.html')  or fnmatch.fnmatch(file,'*.params') or \
           fnmatch.fnmatch(file, '*report.Rmd'):
                os.rename('./'+file,'./report/'+file)
        elif fnmatch.fnmatch(file, '*.pdf'):
                os.rename('./'+file,'./plots/'+file)
        else:
            continue
    except Exception as e:
        print('Error moving/deleting file ' + file + " (" + str(e) + ")")
        errcnt += 1
        continue
if errcnt:
    print("post_analysis encountered " + str(errcnt) + " errors.")
    sys.exit(2)
