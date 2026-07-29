#!/usr/bin/env python3
"""
Validate module.json manifests in WILDS modules and pipelines.

module.json is optional (see .github/CONTRIBUTING.md) and follows the
experimental WDL v1.4 module manifest spec (openwdl/wdl#765). This script
checks structure only where a module.json is present; it does not require
every module/pipeline to have one.
"""

import json
import re
import sys
from pathlib import Path

DEP_NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_-]*$")
TOOL_ID_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9._-]*:.+$")
VERSION_SELECTORS = ("version", "tag", "branch", "commit")


def validate_tool(tool, idx):
    """Validate a single tools[] entry. Returns list of error strings."""
    errors = []
    if not isinstance(tool, dict):
        return [f"  tools[{idx}]: expected an object"]

    for field in ("name", "version", "license"):
        if field not in tool:
            errors.append(f"  tools[{idx}]: missing required field '{field}'")
        elif not isinstance(tool[field], str) or not tool[field]:
            errors.append(f"  tools[{idx}].{field}: expected a non-empty string")

    if "url" in tool and not isinstance(tool["url"], str):
        errors.append(f"  tools[{idx}].url: expected a string")

    if "ids" in tool:
        if not isinstance(tool["ids"], list):
            errors.append(f"  tools[{idx}].ids: expected an array")
        else:
            for cid in tool["ids"]:
                if not isinstance(cid, str) or not TOOL_ID_RE.match(cid):
                    errors.append(
                        f"  tools[{idx}].ids: '{cid}' is not a valid CURIE (expected 'prefix:reference')"
                    )

    for legacy in ("homepage", "doi"):
        if legacy in tool:
            errors.append(
                f"  tools[{idx}]: uses legacy field '{legacy}' — use 'url' for links and "
                f"'ids' (e.g. [\"doi:10.xxxx/...\"]) for identifiers instead"
            )

    return errors


def validate_dependency(name, dep):
    """Validate a single dependencies{} entry. Returns list of error strings."""
    errors = []
    if not DEP_NAME_RE.match(name):
        errors.append(f"  dependencies.{name}: key must match ^[A-Za-z][A-Za-z0-9_-]*$")

    if not isinstance(dep, dict):
        return errors + [f"  dependencies.{name}: expected an object"]

    has_git = "git" in dep
    selectors_present = [s for s in VERSION_SELECTORS if s in dep]

    if has_git:
        if len(selectors_present) != 1:
            errors.append(
                f"  dependencies.{name}: git dependency must specify exactly one of "
                f"'version', 'tag', 'branch', 'commit' (found {selectors_present or 'none'})"
            )
    else:
        if "path" not in dep:
            errors.append(
                f"  dependencies.{name}: must specify either 'git' (with a version selector) or 'path'"
            )
        elif selectors_present:
            errors.append(
                f"  dependencies.{name}: local 'path' dependency must not specify a version selector"
            )

    return errors


def validate_module_json(filepath):
    """Validate a single module.json file. Returns list of error strings."""
    errors = []
    try:
        with open(filepath) as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        return [f"  Invalid JSON: {e}"]

    if not isinstance(data, dict):
        return ["  Expected a JSON object at the top level"]

    for field in ("name", "license"):
        if field not in data:
            errors.append(f"  Missing required field '{field}'")
        elif not isinstance(data[field], str) or not data[field]:
            errors.append(f"  '{field}': expected a non-empty string")

    if "version" in data:
        errors.append(
            "  Unexpected top-level 'version' field — module versioning comes from "
            "Git tags, not the manifest (see openwdl/wdl#765)"
        )

    if "authors" in data and not (
        isinstance(data["authors"], list) and all(isinstance(a, str) for a in data["authors"])
    ):
        errors.append("  'authors': expected an array of strings")

    if "readme" in data and not (data["readme"] is False or isinstance(data["readme"], str)):
        errors.append("  'readme': expected a string path or `false`")

    if "exclude" in data and not (
        isinstance(data["exclude"], list) and all(isinstance(p, str) for p in data["exclude"])
    ):
        errors.append("  'exclude': expected an array of strings")

    module_dir = filepath.parent

    entrypoint = data.get("entrypoint", "index.wdl")
    if isinstance(entrypoint, str) and not (module_dir / entrypoint).exists():
        errors.append(f"  'entrypoint' points to a file that doesn't exist: {entrypoint}")

    readme = data.get("readme", "README.md")
    if isinstance(readme, str) and not (module_dir / readme).exists():
        errors.append(f"  'readme' points to a file that doesn't exist: {readme}")

    if "tools" in data:
        if not isinstance(data["tools"], list):
            errors.append("  'tools': expected an array")
        else:
            for i, tool in enumerate(data["tools"]):
                errors.extend(validate_tool(tool, i))

    if "dependencies" in data:
        if not isinstance(data["dependencies"], dict):
            errors.append("  'dependencies': expected an object")
        else:
            for name, dep in data["dependencies"].items():
                errors.extend(validate_dependency(name, dep))

    return errors


def find_module_json_files(name=None):
    """Find module.json files under modules/ and pipelines/, optionally scoped to one name."""
    files = []
    for tier in ("modules", "pipelines"):
        tier_dir = Path(tier)
        if not tier_dir.is_dir():
            continue
        if name:
            candidate = tier_dir / name / "module.json"
            if candidate.exists():
                files.append(candidate)
        else:
            files.extend(sorted(tier_dir.glob("*/module.json")))
    return files


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Validate module.json manifests in WILDS modules/pipelines.")
    parser.add_argument(
        "name", nargs="?", default=None,
        help="Specific module/pipeline name to validate (e.g., ww-bwa). If omitted, validates all."
    )
    args = parser.parse_args()

    files = find_module_json_files(args.name)

    if not files:
        if args.name:
            print(f"Skipping {args.name} (no module.json found)")
        else:
            print("No module.json files found")
        return 0

    all_errors = {}
    for filepath in files:
        label = f"{filepath.parent.parent.name}/{filepath.parent.name}"
        print(f"Validating {filepath} ...")
        errors = validate_module_json(filepath)
        if errors:
            all_errors[label] = errors
            print(f"  FAIL ({len(errors)} issue(s))")
        else:
            print("  OK")

    if all_errors:
        print(f"\n{'=' * 50}")
        print(f"module.json validation failed for {len(all_errors)} item(s):\n")
        for label, errors in all_errors.items():
            print(f"{label}:")
            for error in errors:
                print(error)
            print()
        return 1

    print(f"\nAll {len(files)} module.json file(s) valid!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
