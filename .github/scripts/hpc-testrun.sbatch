#!/bin/bash
#SBATCH --job-name=wilds_wdl_testrun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/wilds_testrun_%j.out
#SBATCH --error=logs/wilds_testrun_%j.err

# Monthly HPC test run for WILDS WDL Library
# Runs all module/pipeline testruns with sprocket and posts results
# as a comment on a central GitHub tracking issue.
#
# Prerequisites:
#   - gh CLI authenticated (run `gh auth login` once)
#   - sprocket, wdlparse, and Apptainer available via module load
#   - GITHUB_ISSUE_NUMBER set to the tracking issue number
#
# Usage:
#   export GITHUB_ISSUE_NUMBER=<issue_number>
#   export WORK_DIR=/hpc/temp/your-username/wilds-testrun
#   sbatch /path/to/hpc-testrun.sbatch
#
# To create the tracking issue (one-time setup):
#   gh issue create \
#     --repo getwilds/wilds-wdl-library \
#     --title "HPC Monthly Test Run Log" \
#     --body "Rolling log of monthly HPC test runs. Each run posts results as a comment." \
#     --label "testing"

set -eo pipefail

# ---- Configuration ----
REPO="getwilds/wilds-wdl-library"
ISSUE_NUMBER="${GITHUB_ISSUE_NUMBER:?Error: set GITHUB_ISSUE_NUMBER before submitting}"
WORK_DIR="${WORK_DIR:-$(mktemp -d)}"
RUN_DATE="$(date '+%Y-%m-%d %H:%M:%S')"
LOGFILE="$(mktemp)"

# ---- Environment setup ----
# Update these module versions as needed for your HPC
module purge
module load sprocket/0.22.0 Apptainer/1.1.6

# ---- Clone fresh copy of the repo ----
echo "Cloning ${REPO} into ${WORK_DIR}..."
git clone "https://github.com/${REPO}.git" "${WORK_DIR}/wilds-wdl-library"
cd "${WORK_DIR}/wilds-wdl-library"

# ---- Run testruns ----
echo "Starting WILDS WDL test runs at ${RUN_DATE}"
echo "Logging to ${LOGFILE}"

# Run sprocket testruns (continues past failures thanks to Makefile update)
SPROCKET_CONFIG="${SPROCKET_CONFIG:-notes/minimal-sprocket-config.toml}"
if make run_sprocket SPROCKET_CONFIG="${SPROCKET_CONFIG}" 2>&1 | tee "${LOGFILE}"; then
    STATUS="PASSED"
else
    STATUS="FAILED"
fi

# ---- Build GitHub comment ----
# Extract failure summary if any
FAILURES=$(grep -A 100 "=== SPROCKET FAILURES ===" "${LOGFILE}" 2>/dev/null || true)
ITEM_COUNT=$(grep -c "^\.\.\. Running " "${LOGFILE}" 2>/dev/null || echo "0")
FAILED_COUNT=$(grep -c "^\.\.\. FAILED:" "${LOGFILE}" 2>/dev/null || echo "0")
PASSED_COUNT=$((ITEM_COUNT - FAILED_COUNT))

if [ "${STATUS}" = "PASSED" ]; then
    STATUS_EMOJI="white_check_mark"
else
    STATUS_EMOJI="x"
fi

COMMENT_BODY="$(cat <<EOF
## :${STATUS_EMOJI}: HPC Test Run — ${RUN_DATE}

| | |
|---|---|
| **Status** | ${STATUS} |
| **Date** | ${RUN_DATE} |
| **Engine** | sprocket |
| **Items tested** | ${ITEM_COUNT} |
| **Passed** | ${PASSED_COUNT} |
| **Failed** | ${FAILED_COUNT} |
| **SLURM Job ID** | ${SLURM_JOB_ID:-N/A} |
| **Branch** | $(git rev-parse --abbrev-ref HEAD) |
| **Commit** | $(git rev-parse --short HEAD) |

$(if [ -n "${FAILURES}" ]; then
    echo "<details>"
    echo "<summary>Failure Details</summary>"
    echo ""
    echo "\`\`\`"
    echo "${FAILURES}"
    echo "\`\`\`"
    echo "</details>"
fi)
EOF
)"

# ---- Post to GitHub ----
echo ""
echo "Posting results to ${REPO}#${ISSUE_NUMBER}..."
gh issue comment "${ISSUE_NUMBER}" \
    --repo "${REPO}" \
    --body "${COMMENT_BODY}"

echo "Results posted successfully."

# ---- Cleanup ----
rm -f "${LOGFILE}"
rm -rf "${WORK_DIR}/wilds-wdl-library"
