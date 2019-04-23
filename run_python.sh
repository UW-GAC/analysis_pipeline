#!/bin/bash
# run a python script with the following passed in args:
#   1. full path to python script
#   2. python script argument #1
#      ...
#   n. last python argument
f () {
    errcode=$? # save the exit code as the first thing done in the trap function
    echo "run_python: error executing python script $script - error code $errcode"
    exit $errcode  # or use some other value or do return instead
}
trap f ERR
script="$1"
shift 1
echo "Executing $script $@"
$script "$@"
