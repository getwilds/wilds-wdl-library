#!/usr/bin/env python
"""Parse make run_sprocket output and post results to a GitHub tracking issue.

Usage:
    python hpc_testrun_report.py <logfile> <issue_number> [--repo OWNER/REPO]

Reads the log file produced by `make run_sprocket`, builds a markdown summary,
and posts it as a comment on the specified GitHub issue using the `gh` CLI.
"""

import argparse
import os
import re
import subprocess
import sys


def parse_log(logfile: str) -> dict:
    """Parse make run_sprocket output for pass/fail counts and failure details."""
    with open(logfile) as f:
        content = f.read()

    items = re.findall(r"^\.\.\. Running (.+)$", content, re.MULTILINE)
    failures = re.findall(r"^\.\.\. FAILED: (.+)$", content, re.MULTILINE)

    # Extract the failure summary block if present
    failure_block = ""
    match = re.search(r"(=== SPROCKET FAILURES ===.*)", content, re.DOTALL)
    if match:
        failure_block = match.group(1).strip()

    return {
        "total": len(items),
        "failed": len(failures),
        "passed": len(items) - len(failures),
        "failure_names": failures,
        "failure_block": failure_block,
    }


def get_git_info() -> dict:
    """Get current branch and short commit hash."""
    branch = subprocess.check_output(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"], universal_newlines=True
    ).strip()
    commit = subprocess.check_output(
        ["git", "rev-parse", "--short", "HEAD"], universal_newlines=True
    ).strip()
    return {"branch": branch, "commit": commit}


def build_comment(results: dict, run_date: str, slurm_job_id: str, run_type: str = "all") -> str:
    """Build the markdown comment body."""
    git = get_git_info()
    status = "PASSED" if results["failed"] == 0 else "FAILED"
    emoji = "white_check_mark" if status == "PASSED" else "x"
    type_label = run_type if run_type != "all" else "modules + pipelines"

    comment = f"""## :{emoji}: HPC Test Run ({type_label}) — {run_date}

| | |
|---|---|
| **Status** | {status} |
| **Date** | {run_date} |
| **Type** | {type_label} |
| **Engine** | sprocket |
| **Items tested** | {results['total']} |
| **Passed** | {results['passed']} |
| **Failed** | {results['failed']} |
| **SLURM Job ID** | {slurm_job_id} |
| **Branch** | {git['branch']} |
| **Commit** | {git['commit']} |"""

    if results["failure_block"]:
        comment += f"""

<details>
<summary>Failure Details</summary>

```
{results['failure_block']}
```
</details>"""

    return comment


def post_comment(comment: str, issue_number: str, repo: str) -> None:
    """Post the comment to the GitHub issue using the gh REST API."""
    subprocess.run(
        [
            "gh", "api",
            "repos/{}/issues/{}/comments".format(repo, issue_number),
            "-f", "body={}".format(comment),
        ],
        check=True,
    )
    print("Results posted to {}#{}".format(repo, issue_number))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logfile", help="Path to make run_sprocket log output")
    parser.add_argument("issue_number", help="GitHub issue number for the tracking issue")
    parser.add_argument(
        "--repo", default="getwilds/wilds-wdl-library",
        help="GitHub repo in OWNER/REPO format (default: getwilds/wilds-wdl-library)",
    )
    args = parser.parse_args()

    run_date = os.environ.get("RUN_DATE", "unknown")
    slurm_job_id = os.environ.get("SLURM_JOB_ID", "N/A")
    run_type = os.environ.get("TYPE", "all")

    results = parse_log(args.logfile)
    comment = build_comment(results, run_date, slurm_job_id, run_type)
    post_comment(comment, args.issue_number, args.repo)


if __name__ == "__main__":
    main()
