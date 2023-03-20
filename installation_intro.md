# Installation introduction

This document provides an introduction to the pipeline for users who wish to set up a new installation on their own cluster.

## Overview

To run the pipeline, users submit python scripts (e.g., `grm.py`, `null_model.py`, `assoc.py`, etc.).
These python scripts handle submitting a series of R scripts (in the `R/` directory) to a cluster, using a driver bash script (e.g., `runRscript.sh`).
The python scripts also handle job dependencies using the specific cluster's job scheduler.

The pipeline is designed to work with multiple job schedulers and clusters.
Configuration for a specific job schedule is handled by two components:

1. The `TopmedPipeline.py` source code contains a `Cluster` class that can be subclassed to work with a specific job scheduler and cluster. Some pre-defined subclasses are provided, but users can also add new subclasses if one of the existing subclasses does not fit their cluster.
2. The user can provide a cluster configuration file containing job submission options, e.g., the queue to submit to, memory limits, etc. These are stored in JSON format in a file and passed to the python scripts with a keyword argument. An example for the `UW_Cluster` subclass can be found in `cluster_cfg.json`.

## Installation steps

### Initial steps

1. Clone the source code: `git clone https://github.com/UW-GAC/analysis_pipeline.git`

2. Change to the repo directory: `cd analysis_pipeline`

3. Install packages: `R --args "R_library" < install_packages.R`

This installs packages in the `analysis_pipeline/R_library/` directory.

4. If necessary, define a `Cluster` subclass for your cluster. See the Cluster configuration section for more information.


### Cluster configuration in `TopmedPipeline.py`
