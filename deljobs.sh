#!/bin/bash
#    arg1: name of analysis log file
#
if [[ $# -ne 1 ]]; then
    echo "Name of analysis log file is required."
    exit 1
fi
AFILE=$1
# test if slurm or sge
which sbatch > /dev/null 2>&1
if [[ $? -eq 0 ]]; then
    QCMD="squeue -j"
    DCMD="scancel"
else
    which qstat > /dev/null 2>&1
    if [[ $? -eq 0 ]]; then
        QCMD="qstat -j"
        DCMD="qdel"
    else
        echo "deljobs.sh only supported on SGE or slurm"
        exit 1
    fi
fi
# get the jobids
JIDS=$(grep "jobid:" $AFILE | awk {'print $NF'})

# iterate and delete the jobs that exist
for JOB in $JIDS
do
    $QCMD $JOB > /dev/null 2>&1
    if [[ $? -eq 0 ]]; then
        $DCMD $JOB > /dev/null 2>&1
        if [[ $? -eq 0 ]]; then
            echo "Deleted job $JOB"
        fi
    else
        echo "Job $JOB doesn't exist"
    fi
done
