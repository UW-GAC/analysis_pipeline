#! /bin/bash

# arguments: sample_list in_file_prefix in_file_suffix out_file_prefix out_file_suffix

if [ "$SGE_TASK_ID" == "23" ]; then
    CHR="X"
else
    CHR=$SGE_TASK_ID
fi

in_file=$2${CHR}$3
out_file=$4${CHR}$5

bcftools view --samples-file $1 -O z -o $out_file $in_file
