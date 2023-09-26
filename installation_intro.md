# Installation introduction

This document provides an introduction to the pipeline for users who wish to set up a new installation on their own cluster.

## Overview

To run the pipeline, users submit python scripts (e.g., `grm.py`, `null_model.py`, `assoc.py`, etc.).
These python scripts handle submitting a series of R scripts (in the `R/` directory) to a cluster, using a driver bash script (e.g., `runRscript.sh`).
The python scripts also handle job dependencies using the specific cluster's job scheduler.

The pipeline is designed to work with multiple job schedulers and clusters.
Configuration for a specific job schedule is handled by two components:

1. The `TopmedPipeline.py` source code contains a `Cluster` class that can be subclassed to work with a specific job scheduler and cluster. Some pre-defined subclasses are provided, but users can also add new subclasses if one of the existing subclasses does not fit their cluster.
2. The user can provide a cluster configuration file containing job submission options, e.g., the queue to submit to, memory limits, etc. These are stored in JSON format in a file and passed to the python scripts with a keyword argument. An example for the `GAC_Slurm_Cluster` subclass can be found in `slurm_cluster_cfg.json`.

## Installation steps

### Python and R installation

You will need a recent installation of R and an installation of python3.

We recommend building R with [Intel MKL](https://software.intel.com/en-us/intel-mkl) for improved performance in PC-Relate and association tests.

### Initial steps

1. Clone the source code: `git clone https://github.com/UW-GAC/analysis_pipeline.git`

3. Install required packages: `R --args "<path_to_R_library>" < install_packages.R`.

4. Install other dependencies as needed. See the "Install dependencies" section for more information.

5. If necessary, define a `Cluster` subclass for your cluster. See the Cluster configuration section for more information.

### Install external dependencies

Some of the scripts have dependencies on external programs. You may need to install the following:

- [bcftools](http://www.htslib.org/download/) - used by `vcf2gds.py`
- [PLINK](https://www.cog-genomics.org/plink2/) - used by `king.py`
- [KING 2.2.4](https://www.kingrelatedness.com/executables/Linux-king224.tar.gz) - used by `king.py`
- [LocusZoom](https://github.com/UW-GAC/locuszoom-standalone) - used by `locuszoom.py`
- [pandoc](https://pandoc.org/installing.html) - required for `null_model.py` and `assoc.py` reports

These will need to be on your path if you want to use the python script that requires them.


### Cluster configuration in `TopmedPipeline.py`

If there is not an appropriate `Cluster` subclass that works with your current job scheduler, you will need to define a new one.

You will need to create a subclass of the `Cluster` class in `TopmedPipeline.py`.
Your subclass will typically need to override the `__init__`, `submitJob`, and `runCommand` methods.

#### The `__init__` method

The `__init__` method needs to set the `class_name` and `std_cluster_file` attributes.
Other tasks can be performed here if necessary for general cluster setup.

#### The `runCommand` method

This method handles running a command without submitting it to the cluster.
It should set the environment variables (e.g., `R_LIBS`, `PATH`) properly.

#### The `submitJob` method

This method is intended to handle converting the option specified in your cluster configuration file to specific arguments to your job scheduler.

For example, in SGE, you submit a job and specify the name with the command `qsub -N <job_name> ...`.
In Slurm, the same command is `sbatch --job-name=<job_name>`.
The `submitJob` method is responsible for processing options specified in the python scripts (e.g., `null_model.py`) and converting them to the appropriate arguments for your job scheduler.

The `submitJob` method should plan to handle the following **kwargs:

- job_name
- hold_array or holdid
- array_range
- request_cores
- memory_limits
- email
- args
- print_only
- enable_resume
- XXX others?

### Cluster configuration JSON file

You can modify the `slurm_cluster_cfg.json` or `sge_cluster_cfg.json` files as needed for your cluster.

Here is a brief description of the various keys:

- `"memory_limits"`: specify the amount of memory that should be requested for each job.
- `"submit_cmd"`: specify the command used to submit a job to the cluster (e.g., `sbatch` for Slurm or `qsub` for SGE).
- `"submit_opts"`: specify any additional options that should be passed to the `submit_cmd`, e.g., specifying the partition or queue, email parameters, etc. Note that the python scripts and the `Cluster` subclass you are using also set parameters when submitting a command.

## Some tips and troubleshooting

- Make sure that you set your R library path to include the packages installed by `install_packages.py`.
- Make sure that your path (e.g., the `$PATH` environment variable) includes the programs you installed in the "Install external dependencies" section.
