#! /bin/bash

#$ -cwd
#$ -S /bin/bash
#$ -j y

args=("$@") # all arguments
unset args[0] # remove first argument (R script name)

export R_LIBS=/projects/geneva/gcc-fs2/R_packages/library

if [ "$SGE_TASK_ID" == "undefined" ] || [ "$SGE_TASK_ID" == "" ]; then
    CHR=""
else 
    CHR="--chromosome $SGE_TASK_ID"
fi

R -q --vanilla --args ${args[@]} $CHR < $1
