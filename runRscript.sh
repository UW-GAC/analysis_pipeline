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

export R_LIBS=/projects/resources/gactools/R_packages/library

R-3.3.2 -q --vanilla --args ${args[@]} $TASK < $1
