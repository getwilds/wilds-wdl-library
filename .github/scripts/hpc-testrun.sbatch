#!/bin/bash
#SBATCH --job-name=wilds_wdl_testrun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --output=logs/wilds_testrun_%j.out
#SBATCH --error=logs/wilds_testrun_%j.err

# Monthly HPC test run for WILDS WDL Library
# Runs all module/pipeline testruns with sprocket and posts results
# as a comment on a central GitHub tracking issue.
#
# Prerequisites:
#   - gh CLI authenticated (run `gh auth login` once)
#   - sprocket and Apptainer available via module load
#   - GITHUB_ISSUE_NUMBER set to the tracking issue number
#
# Usage:
#   export GITHUB_ISSUE_NUMBER=<issue_number>
#   export WORK_DIR=/path/to/your/working/directory
#   sbatch /path/to/hpc-testrun.sbatch
#
# Optional:
#   BRANCH=my-feature-branch  # Defaults to "main"
#   TYPE=modules|pipelines    # Defaults to "all"
#
# To create the tracking issue (one-time setup):
#   gh issue create \
#     --repo getwilds/wilds-wdl-library \
#     --title "HPC Monthly Test Run Log" \
#     --body "Rolling log of monthly HPC test runs. Each run posts results as a comment." \
#     --label "testing"

set -eo pipefail

# ---- Configuration ----
ISSUE_NUMBER="${GITHUB_ISSUE_NUMBER:?Error: set GITHUB_ISSUE_NUMBER before submitting}"
WORK_DIR="${WORK_DIR:-$(mktemp -d)}"
export RUN_DATE="$(date '+%Y-%m-%d %H:%M:%S')"
LOGFILE="$(mktemp)"
REPO="getwilds/wilds-wdl-library"

# ---- Environment setup ----
# Update these module versions as needed for your HPC
module purge
module load sprocket/0.23.0 Apptainer/1.1.6

# ---- Clone fresh copy of the repo ----
BRANCH="${BRANCH:-main}"
CLONE_DIR="${WORK_DIR}/wilds-wdl-library-${SLURM_JOB_ID}"
echo "Cloning ${REPO} (branch: ${BRANCH}) into ${CLONE_DIR}..."
git clone --branch "${BRANCH}" "https://github.com/${REPO}.git" "${CLONE_DIR}"
cd "${CLONE_DIR}"

# ---- Run testruns ----
export TYPE="${TYPE:-all}"
echo "Starting WILDS WDL test runs at ${RUN_DATE} (type: ${TYPE})"
SPROCKET_CONFIG="${SPROCKET_CONFIG:-.github/configs/sprocket-hpc.toml}"
make run_sprocket TYPE="${TYPE}" SPROCKET_CONFIG="${SPROCKET_CONFIG}" SPROCKET_VERBOSITY="-q" 2>&1 | tee "${LOGFILE}" || true

# ---- Post results to GitHub ----
python3 .github/scripts/hpc_testrun_report.py "${LOGFILE}" "${ISSUE_NUMBER}" --repo "${REPO}"

# ---- Cleanup ----
rm -f "${LOGFILE}"
# Note: cloned repo + sprocket output left in place for debugging
