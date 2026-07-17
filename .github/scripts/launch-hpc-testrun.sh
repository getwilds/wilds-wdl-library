#!/bin/bash
# Launch script for monthly HPC test runs.
# Intended to be called from crontab or manually.
#
# Update the paths and issue number below for your environment.
#
# Crontab example (1st of every month at midnight):
#   0 0 1 * * /path/to/launch-hpc-testrun.sh

export GITHUB_ISSUE_NUMBER=324
export WORK_DIR=/path/to/scratch/wilds-testrun

SBATCH_SCRIPT=/path/to/hpc-testrun.sbatch

# Submit modules and pipelines as separate SLURM jobs
TYPE=modules NUM_RETRIES=2 sbatch "${SBATCH_SCRIPT}"
TYPE=pipelines NUM_RETRIES=2 sbatch "${SBATCH_SCRIPT}"
