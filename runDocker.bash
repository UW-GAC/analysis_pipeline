#!/bin/bash
# a bash script invoked by sbatch in slurm to insure topmed's nfs data volume is
# mounted and then executes the python script to execute an analysis pipeline
# command.
#
# arguments:
#   1 python submit script (full path)
#   2-n arguments passed to the python script
f () {
    errcode=$? # save the exit code as the first thing done in the trap function
    echo "runDocker.bash: error executing python script $PSCRIPT - error code $errcode"
    exit $errcode  # or use some other value or do return instead
}
trap f ERR
PSCRIPT=$1
shift 1
CDIR=/projects/topmed
STIME=10
MAXTIME=180
CTR=0
# loop to check for mounted directory
while :
do
    if [ ! -d $CDIR ]; then
        sleep $STIME
        ((CTR=CTR+STIME))
        if [ $CTR -ge $MAXTIME ]; then
            echo "Error: check directory $CDIR could not be found."
            exit 2
        fi
    else
        break
    fi
done
shift 1
echo "Check directory $CDIR exists; running python script:"
echo "> $PSCRIPT $@"
#
$PSCRIPT "$@"
