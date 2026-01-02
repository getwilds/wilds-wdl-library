#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   discover_wdls.py
@Time    :   2025/10/16
@Author  :   Taylor Firman
@Version :   v1.0
@Contact :   tfirman@fredhutch.org
@Desc    :   Unified script to identify WILDS modules or pipelines needing test runs based on PR changes.
'''

import json
import os
import subprocess
import sys
from pathlib import Path

# Valid item types
VALID_TYPES = ['modules', 'pipelines']

def get_changed_files():
    """Get list of changed files in PR, or None for workflow_dispatch"""
    if os.environ.get('GITHUB_EVENT_NAME') == 'pull_request':
        try:
            result = subprocess.run([
                'git', 'diff', '--name-only',
                f"origin/{os.environ['GITHUB_BASE_REF']}...HEAD"
            ], capture_output=True, text=True, check=True)
            return result.stdout.strip().split('\n') if result.stdout.strip() else []
        except subprocess.CalledProcessError as e:
            print(f"ERROR: Could not get changed files: {e}")
            print("This is a critical error - we cannot safely determine what to test")
            raise SystemExit(1)
    return None

def find_valid_items(item_type):
    """Find all directories with valid structure (containing testrun.wdl)"""
    items = []
    items_dir = Path(item_type)
    singular = item_type.rstrip('s')  # Remove trailing 's' for singular

    if not items_dir.exists():
        print(f"No {item_type} directory found")
        return items

    for item_dir in items_dir.iterdir():
        if not item_dir.is_dir():
            continue

        item_name = item_dir.name
        wdl_file = item_dir / 'testrun.wdl'

        if wdl_file.exists():
            items.append(item_name)
            print(f"Found valid {singular}: {item_name}")
        else:
            print(f"Skipping {item_name} - missing testrun.wdl")

    return sorted(items)

def filter_items_by_changes(items, changed_files, item_type):
    """Filter items to only those with WDL file changes"""
    if changed_files is None:
        print("No file filtering applied (potential git error), skipping...")
        return []

    if not changed_files:
        print("No files changed in this PR")
        return []

    print(f"Files changed in PR: {len(changed_files)}")
    for f in changed_files[:10]:  # Show first 10
        print(f"  {f}")
    if len(changed_files) > 10:
        print(f"  ... and {len(changed_files) - 10} more")

    filtered = []
    for item in items:
        # Only consider .wdl file changes
        item_wdl_files = [f for f in changed_files
                          if f.startswith(f'{item_type}/{item}/') and f.endswith('.wdl')]
        if item_wdl_files:
            filtered.append(item)
            print(f"Including {item} ({len(item_wdl_files)} changed .wdl files)")
        else:
            print(f"Skipping {item} (no .wdl file changes)")

    return filtered

def handle_workflow_dispatch_input(items, dispatch_item, item_type):
    """Handle manual workflow dispatch input for specific item"""
    if dispatch_item is None:
        dispatch_item = ''
    dispatch_item = dispatch_item.strip()
    singular = item_type.rstrip('s')  # Remove trailing 's' for singular

    if not dispatch_item:
        print(f"No specific {singular} requested for workflow_dispatch, running all {item_type}")
        return items

    # Check if the requested item exists and is valid
    if dispatch_item in items:
        print(f"Running specific {singular} from workflow_dispatch: {dispatch_item}")
        return [dispatch_item]
    else:
        print(f"ERROR: Requested {singular} '{dispatch_item}' not found or invalid")
        print(f"Available {item_type}: {', '.join(items)}")
        # Return empty list to fail gracefully rather than running wrong items
        return []

def main():
    # Parse command line arguments
    if len(sys.argv) < 2:
        print("ERROR: Item type not specified")
        print("Usage: discover_wdls.py <modules|pipelines> [specific_item_name]")
        raise SystemExit(1)

    item_type = sys.argv[1].lower()
    if item_type not in VALID_TYPES:
        print(f"ERROR: Invalid item type '{item_type}'")
        print(f"Valid types: {', '.join(VALID_TYPES)}")
        raise SystemExit(1)

    singular = item_type.rstrip('s')  # Remove trailing 's' for singular

    dispatch_item = None
    if len(sys.argv) > 2:
        dispatch_item = sys.argv[2]
        print(f"Workflow dispatch {singular} argument: '{dispatch_item}'")

    print(f"=== WILDS {item_type.title()} Discovery ===")

    # Find all valid items
    all_items = find_valid_items(item_type)
    print(f"\nFound {len(all_items)} valid {item_type}")

    # Check event type
    event_name = os.environ.get('GITHUB_EVENT_NAME', '')
    print(f"GitHub event: {event_name}")

    if event_name == 'workflow_dispatch':
        # Handle manual workflow dispatch
        final_items = handle_workflow_dispatch_input(all_items, dispatch_item, item_type)
    else:
        # Handle PR - filter by changes
        changed_files = get_changed_files()
        final_items = filter_items_by_changes(all_items, changed_files, item_type)

    print(f"\nFinal selection: {len(final_items)} {item_type}")
    for item in final_items:
        print(f"  - {item}")

    # Output for GitHub Actions
    items_json = json.dumps(final_items)

    # Write to GITHUB_OUTPUT using item_type as the key
    github_output = os.environ.get('GITHUB_OUTPUT')
    if github_output:
        with open(github_output, 'a') as f:
            f.write(f"{item_type}={items_json}\n")
    else:
        print(f"{item_type}={items_json}")

if __name__ == "__main__":
    main()
