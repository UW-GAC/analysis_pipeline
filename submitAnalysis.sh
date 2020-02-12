#!/bin/bash
# a bash script invoked by sbatch in slurm as a preprocess step before executing
# the python script for executing the pipeline.
#
# arguments:
#   1 function name
#   2-n function specific args (0 to n)
#   n+1 python script (full path)
#   n+2-m python script args
#
# check mountpoint args:
#   none
# create label args:
#   1 label
#   2 label value
# none:
#   none
err_trap () {
    errcode=$? # save the exit code as the first thing done in the trap function
    echo "PrePython>>> Error on line $1 - error code $errcode"
    exit $errcode  # or use some other value or do return instead
}
trap 'err_trap $LINENO' ERR
# check mount point
check_mountpoint() {
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
                echo "PrePython>>> Error: check directory $CDIR could not be found."
                exit 2
            fi
        else
            break
        fi
    done
    echo "PrePython> Check directory $CDIR exists; running python script:"
    SHIFTS=0
}
create_label() {
    # get the label and shift
    LABEL=$1
    VALUE=$2
    # get the instance name
    INAME=`curl  http://metadata.google.internal/computeMetadata/v1/instance/name -H "Metadata-Flavor: Google"`
    # create the label on the instance
    echo "PrePython> Creating label $LABEL on instance $INAME"
    gcloud compute instances add-labels $INAME --labels=$LABEL=$VALUE --zone us-west1-a
    SHIFTS=2
}
none() {
    SHIFTS=0
}
FUNC=$1
shift 1
$FUNC $@
echo "PrePython> Shift $SHIFTS to get to python script"
shift $SHIFTS
PSCRIPT=$1
shift 1
echo "PrePython> $PSCRIPT $@"
#
$PSCRIPT "$@"
