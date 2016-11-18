#! /bin/bash

#$ -cwd
#$ -S /bin/bash
#$ -j y

while getopts ":cs" opt; do
  case $opt in
    c)
      TASK="--chromosome $SGE_TASK_ID"
      ;;
    s)
      TASK="--segment $SGE_TASK_ID"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done
shift "$((OPTIND-1))"

if [ "$SGE_TASK_ID" == "undefined" ] || [ "$SGE_TASK_ID" == "" ]; then
    TASK=""
fi

args=("$@") # all arguments
unset args[0] # remove first argument (R script name)

export R_LIBS=/projects/geneva/gcc-fs2/R_packages/library
#temporary for testing
#export R_LIBS=/projects/users/stephanie/Code/TOPMed/analysis_pipeline/R_library:/projects/geneva/gcc-fs2/R_packages/library

R -q --vanilla --args ${args[@]} $TASK < $1
