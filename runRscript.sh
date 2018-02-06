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
tb=`date`
echo ">>> Start job: $JOB_ID at $tb ... "
R -q --vanilla --args ${args[@]} $TASK < $1


export R_exit_code="$?"

if [ $R_exit_code -ne "0" ]
then
  if [[ "$SGE_TASK_ID" -ne "" ]]; then
    touch fail.${JOB_ID}.${SGE_TASK_ID}
  else
    touch fail.${JOB_ID}
  fi
fi
te=`date`
echo ">>> End job: $JOB_ID at $te"
exit $R_exit_code
