#! /bin/bash

# TODO: UPDATE FOR SLURM?
# Set cross-cluster environment variables.


# Determine cluster type and set cluster-agnostic variables appropriate.
if [[ ! -z "${SGE_ROOT}" ]]
then
  CLUSTER_JOB_ID=$JOB_ID
  CLUSTER_TASK_ID=$SGE_TASK_ID
elif [[ ! -z "${SLURM_CLUSTER_NAME}" ]]
then
  export CLUSTER_JOB_ID=$SLURM_JOB_ID
  export CLUSTER_TASK_ID=$SLURM_ARRAY_TASK_ID
fi

# set MKL_NUM_THREADS to match number of available cores
export MKL_NUM_THREADS=$CLUSTER_SLOTS

while getopts ":cs" opt; do
  case $opt in
    c)
      TASK="--chromosome $CLUSTER_TASK_ID"
      ;;
    s)
      TASK="--segment $CLUSTER_TASK_ID"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done
shift "$((OPTIND-1))"

if [ "$CLUSTER_TASK_ID" == "undefined" ] || [ "$CLUSTER_TASK_ID" == "" ]; then
    TASK=""
fi

args=("$@") # all arguments
unset args[0] # remove first argument (R script name)
tb=`date`
echo ">>> Start job: $CLUSTER_JOB_ID at $tb ... "
R -q --vanilla --args ${args[@]} $TASK < $1


R_exit_code="$?"
r_status=$R_exit_code
if [[ $R_exit_code -ne "0" ]]; then
  echo ">>> Error: R status code $R_exit_code"
  if [[ "$CLUSTER_TASK_ID" -ne "" ]]; then
    touch fail.${CLUSTER_JOB_ID}.${CLUSTER_TASK_ID}
  else
    touch fail.${CLUSTER_JOB_ID}
  fi
  if [[ ! -z "${SGE_ROOT+x}" && ! -z "${ENABLE_EQW+x}" ]]; then
    r_status=100
    echo ">>> Error: return status $r_status"
  fi
fi

te=`date`
echo ">>> End job: $CLUSTER_JOB_ID at $te"
exit $r_status
