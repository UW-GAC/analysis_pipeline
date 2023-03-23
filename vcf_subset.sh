#! /bin/bash

# arguments: sample_list in_file_prefix in_file_suffix out_file_prefix out_file_suffix

which sbatch > /dev/null 2>&1
if [[ $? -eq 0 ]]; then
  export CLUSTER_JOB_ID=$SLURM_JOB_ID
  export CLUSTER_TASK_ID=$SLURM_ARRAY_TASK_ID
  export NSLOTS=$SLURM_JOB_CPUS_PER_NODE
else
  which qstat > /dev/null 2>&1
  if [[ $? -eq 0 ]]; then
    CLUSTER_JOB_ID=$JOB_ID
    CLUSTER_TASK_ID=$SGE_TASK_ID
  else
      echo "runRscript.sh only supported on SGE or slurm"
      exit 1
  fi
fi


if [ "$CLUSTER_TASK_ID" == "23" ]; then
    CHR="X"
else
    CHR=$CLUSTER_TASK_ID
fi

in_file=$2${CHR}$3
out_file=$4${CHR}$5

bcftools view --min-ac 1 --samples-file $1 -O z -o $out_file $in_file
