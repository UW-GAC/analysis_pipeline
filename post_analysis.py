#! /usr/bin/env python
import     time
from       argparse import ArgumentParser
import     datetime
import     fnmatch
import     os

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
for file in os.listdir('.'):
    try:
        if fnmatch.fnmatch(file, '*.log') or fnmatch.fnmatch(file,'*.trace') or \
           fnmatch.fnmatch(file, '*.o*') or fnmatch.fnmatch(file,'*.po*'):
                os.rename('./'+file,'./log/'+file)
        if fnmatch.fnmatch(file, '*report.*')  or fnmatch.fnmatch(file,'*.params'):
            os.rename('./'+file,'./report/'+file)
        if fnmatch.fnmatch(file, '*.pdf'):
            os.rename('./'+file,'./plots/'+file)
    except Exception as e:
        print('Error moving file ' + file + " (" + str(e) + ")")
        errcnt += 1
        continue
if errcnt:
    print("post_analysis encountered " + str(errcnt) + " errors.")
    sys.exit(2)
