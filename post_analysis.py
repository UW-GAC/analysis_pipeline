#! /usr/bin/env python
import     time
from       argparse import ArgumentParser
from       datetime import datetime, timedelta
import     fnmatch
import     os


# parse input
parser = ArgumentParser( description = "script for post analysis stuff" )
parser.add_argument( "-a", "--analysis",
                     help = "Name of analysis" )
parser.add_argument( "-l", "--logfile",
                     help = "Name of analysis log file" )
parser.add_argument( "-s", "--starttime",
                     help = "start time (in seconds) of analysis ")

args = parser.parse_args()
analysis = args.analysis
starttime = float(args.starttime)
logfile = args.logfile

end_str = time.asctime()
tFmt = '%a %b %d %H:%M:%S %Y'
dt_ref = datetime(1970,1,1)
end_time = (datetime.strptime(end_str,tFmt) - dt_ref).total_seconds()
td_hrs = (end_time-starttime)/60./60.
# append to log file
with open(logfile, "a") as afile:
    afile.write(analysis + " end time: " + end_str + "\n")
    afile.write(analysis + " elapse time(hrs): " + str(td_hrs) + "\n")

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
