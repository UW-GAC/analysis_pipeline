#! /bin/bash

# set MKL_NUM_THREADS to match number of available cores
export MKL_NUM_THREADS=$NSLOTS

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


R_exit_code="$?"
r_status=$R_exit_code
if [[ $R_exit_code -ne "0" ]]; then
  echo ">>> Error: R status code $R_exit_code"
  if [[ "$SGE_TASK_ID" -ne "" ]]; then
    touch fail.${JOB_ID}.${SGE_TASK_ID}
  else
    touch fail.${JOB_ID}
  fi
  if [[ ! -z "${SGE_ROOT+x}" && ! -z "${ENABLE_EQW+x}" ]]; then
    r_status=100
    echo ">>> Error: return status $r_status"
  fi
fi

te=`date`
echo ">>> End job: $JOB_ID at $te"
exit $r_status
