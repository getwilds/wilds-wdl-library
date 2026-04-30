"""Collect GitHub repository traffic stats and merge them into CSV history.

Reads the views/clones/referrers/paths endpoints for the repo named in the
GITHUB_REPOSITORY env var and writes results into ./data/. Daily-granularity
endpoints (views, clones) are merged into a single CSV per metric, with
later snapshots overwriting earlier rows for the same date. Top-10
endpoints (referrers, paths) only return a current snapshot, so each run
writes a dated JSON file under data/referrers/ and data/paths/.

Auth uses the GH_TOKEN env var (set to GITHUB_TOKEN in the workflow).
"""

import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd

API_ROOT = "https://api.github.com"
DATA_DIR = Path("data")


def gh_get(path: str) -> dict:
    token = os.environ["GH_TOKEN"]
    req = Request(
        f"{API_ROOT}/{path}",
        headers={
            "Authorization": f"Bearer {token}",
            "Accept": "application/vnd.github+json",
            "X-GitHub-Api-Version": "2022-11-28",
        },
    )
    with urlopen(req) as resp:
        return json.load(resp)


def merge_daily(metric: str, repo: str) -> None:
    """Fetch a daily-granularity metric and merge it into data/<metric>.csv."""
    payload = gh_get(f"repos/{repo}/traffic/{metric}")
    fresh = pd.DataFrame(payload[metric])
    if not fresh.empty:
        fresh["date"] = fresh["timestamp"].str[:10]
        fresh = fresh[["date", "count", "uniques"]]

    csv_path = DATA_DIR / f"{metric}.csv"
    if csv_path.exists():
        history = pd.read_csv(csv_path, dtype={"date": str})
        # Drop overlapping dates from history so the fresh snapshot wins
        history = history[~history["date"].isin(fresh["date"])]
        merged = pd.concat([history, fresh], ignore_index=True)
    else:
        merged = fresh

    merged = merged.sort_values("date").reset_index(drop=True)
    merged.to_csv(csv_path, index=False)


def snapshot_top(metric: str, repo: str, today: str) -> None:
    """Fetch a top-N metric and write a dated JSON snapshot."""
    out_dir = DATA_DIR / metric
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = gh_get(f"repos/{repo}/traffic/popular/{metric}")
    (out_dir / f"{today}.json").write_text(json.dumps(payload, indent=2) + "\n")


def plot_metric(metric: str) -> None:
    """Render a simple line chart of count + uniques over time."""
    csv_path = DATA_DIR / f"{metric}.csv"
    if not csv_path.exists():
        return
    df = pd.read_csv(csv_path, parse_dates=["date"])
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(df["date"], df["count"], label="Total", marker="o", markersize=3)
    ax.plot(df["date"], df["uniques"], label="Unique", marker="o", markersize=3)
    ax.set_title(f"{metric.capitalize()} over time")
    ax.set_ylabel(metric.capitalize())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    fig.autofmt_xdate()
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(DATA_DIR / f"{metric}.png", dpi=120)
    plt.close(fig)


def main() -> int:
    repo = os.environ["GITHUB_REPOSITORY"]
    today = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    DATA_DIR.mkdir(exist_ok=True)

    merge_daily("views", repo)
    merge_daily("clones", repo)
    snapshot_top("referrers", repo, today)
    snapshot_top("paths", repo, today)

    plot_metric("views")
    plot_metric("clones")
    return 0


if __name__ == "__main__":
    sys.exit(main())
