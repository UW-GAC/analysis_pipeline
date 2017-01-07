#! /bin/bash

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

R -q --vanilla --args ${args[@]} $TASK < $1
