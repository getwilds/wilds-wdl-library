#!/usr/bin/env python3
"""
Sign module.json manifests in WILDS modules and pipelines, writing/updating
each one's module.sig.

Only re-signs a module/pipeline whose directory content actually changed
between two given git refs, so CI doesn't spend time re-invoking `sprocket
dev module sign` on every module on every push. Pass --all to sign every
module/pipeline with a module.json regardless of what changed.

Requires a `sprocket` build with `dev module sign` (not yet in a released
Sprocket version as of this writing; see .github/workflows/sign-modules.yml
for which branch/commit CI installs).
"""

import argparse
import subprocess
import sys
from pathlib import Path


def find_manifest_dirs():
    """Find every modules/*/ and pipelines/*/ directory containing a module.json."""
    dirs = []
    for tier in ("modules", "pipelines"):
        tier_dir = Path(tier)
        if not tier_dir.is_dir():
            continue
        for item_dir in sorted(tier_dir.iterdir()):
            if (item_dir / "module.json").exists():
                dirs.append(item_dir)
    return dirs


def changed_manifest_dirs(before_ref, after_ref):
    """Return manifest dirs whose contents differ between before_ref and after_ref."""
    result = subprocess.run(
        ["git", "diff", "--name-only", before_ref, after_ref],
        capture_output=True, text=True, check=True,
    )
    changed_top_dirs = set()
    for line in result.stdout.splitlines():
        parts = line.split("/", 2)
        if len(parts) >= 2 and parts[0] in ("modules", "pipelines"):
            changed_top_dirs.add(Path(parts[0]) / parts[1])

    return sorted(d for d in changed_top_dirs if (d / "module.json").exists())


def sign_manifest_dir(manifest_dir, key_path):
    """Run `sprocket dev module sign` against a single manifest directory."""
    subprocess.run(
        [
            "sprocket", "dev", "module", "sign",
            "--key", str(key_path),
            "--manifest-path", str(manifest_dir),
        ],
        check=True,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Sign module.json manifests in WILDS modules/pipelines."
    )
    parser.add_argument("--key", required=True, help="Path to the OpenSSH-format Ed25519 private key")
    parser.add_argument(
        "--all", action="store_true",
        help="Sign every module/pipeline, not just those changed between --before and --after"
    )
    parser.add_argument("--before", default=None, help="Git ref to diff from (ignored with --all)")
    parser.add_argument("--after", default="HEAD", help="Git ref to diff to (ignored with --all; default: HEAD)")
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print which modules/pipelines would be signed without actually signing them"
    )
    args = parser.parse_args()

    if args.all:
        targets = find_manifest_dirs()
        print(f"Signing all {len(targets)} module(s)/pipeline(s) with a module.json")
    else:
        if not args.before:
            parser.error("--before is required unless --all is passed")
        targets = changed_manifest_dirs(args.before, args.after)
        print(f"Signing {len(targets)} module(s)/pipeline(s) changed between {args.before} and {args.after}")

    if not targets:
        print("Nothing to sign")
        return 0

    for target in targets:
        if args.dry_run:
            print(f"  [dry run] would sign {target}")
            continue
        print(f"  Signing {target} ...")
        sign_manifest_dir(target, args.key)

    return 0


if __name__ == "__main__":
    sys.exit(main())
