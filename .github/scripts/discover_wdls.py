#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   discover_wdls.py
@Time    :   2025/10/16
@Author  :   Taylor Firman
@Version :   v1.0
@Contact :   tfirman@fredhutch.org
@Desc    :   Unified script to identify WILDS modules or vignettes needing test runs based on PR changes.
'''

import json
import os
import subprocess
import sys
from pathlib import Path

# Configuration for different item types
ITEM_CONFIGS = {
    'modules': {
        'directory': 'modules',
        'wdl_pattern': lambda name: f"{name}.wdl",  # modules/{name}/{name}.wdl
        'output_key': 'modules',
        'singular': 'module',
        'plural': 'modules'
    },
    'vignettes': {
        'directory': 'vignettes',
        'wdl_pattern': lambda name: 'testrun.wdl',  # vignettes/{name}/testrun.wdl
        'output_key': 'vignettes',
        'singular': 'vignette',
        'plural': 'vignettes'
    }
}

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

def find_valid_items(config):
    """Find all directories with valid structure based on config"""
    items = []
    items_dir = Path(config['directory'])

    if not items_dir.exists():
        print(f"No {config['directory']} directory found")
        return items

    for item_dir in items_dir.iterdir():
        if not item_dir.is_dir():
            continue

        item_name = item_dir.name
        wdl_file = item_dir / config['wdl_pattern'](item_name)

        if wdl_file.exists():
            items.append(item_name)
            print(f"Found valid {config['singular']}: {item_name}")
        else:
            print(f"Skipping {item_name} - missing {wdl_file.name}")

    return sorted(items)

def filter_items_by_changes(items, changed_files, config):
    """Filter items to only those with changes"""
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
        item_files = [f for f in changed_files if f.startswith(f'{config["directory"]}/{item}/')]
        if item_files:
            filtered.append(item)
            print(f"Including {item} ({len(item_files)} changed files)")
        else:
            print(f"Skipping {item} (no changes)")

    return filtered

def handle_workflow_dispatch_input(items, dispatch_item, config):
    """Handle manual workflow dispatch input for specific item"""
    if dispatch_item is None:
        dispatch_item = ''
    dispatch_item = dispatch_item.strip()

    if not dispatch_item:
        print(f"No specific {config['singular']} requested for workflow_dispatch, running all {config['plural']}")
        return items

    # Check if the requested item exists and is valid
    if dispatch_item in items:
        print(f"Running specific {config['singular']} from workflow_dispatch: {dispatch_item}")
        return [dispatch_item]
    else:
        print(f"ERROR: Requested {config['singular']} '{dispatch_item}' not found or invalid")
        print(f"Available {config['plural']}: {', '.join(items)}")
        # Return empty list to fail gracefully rather than running wrong items
        return []

def main():
    # Parse command line arguments
    if len(sys.argv) < 2:
        print("ERROR: Item type not specified")
        print("Usage: discover_wdls.py <modules|vignettes> [specific_item_name]")
        raise SystemExit(1)

    item_type = sys.argv[1].lower()
    if item_type not in ITEM_CONFIGS:
        print(f"ERROR: Invalid item type '{item_type}'")
        print(f"Valid types: {', '.join(ITEM_CONFIGS.keys())}")
        raise SystemExit(1)

    config = ITEM_CONFIGS[item_type]

    dispatch_item = None
    if len(sys.argv) > 2:
        dispatch_item = sys.argv[2]
        print(f"Workflow dispatch {config['singular']} argument: '{dispatch_item}'")

    print(f"=== WILDS {config['plural'].title()} Discovery ===")

    # Find all valid items
    all_items = find_valid_items(config)
    print(f"\nFound {len(all_items)} valid {config['plural']}")

    # Check event type
    event_name = os.environ.get('GITHUB_EVENT_NAME', '')
    print(f"GitHub event: {event_name}")

    if event_name == 'workflow_dispatch':
        # Handle manual workflow dispatch
        final_items = handle_workflow_dispatch_input(all_items, dispatch_item, config)
    else:
        # Handle PR - filter by changes
        changed_files = get_changed_files()
        final_items = filter_items_by_changes(all_items, changed_files, config)

    print(f"\nFinal selection: {len(final_items)} {config['plural']}")
    for item in final_items:
        print(f"  - {item}")

    # Output for GitHub Actions
    items_json = json.dumps(final_items)

    # Write to GITHUB_OUTPUT using the appropriate key
    github_output = os.environ.get('GITHUB_OUTPUT')
    if github_output:
        with open(github_output, 'a') as f:
            f.write(f"{config['output_key']}={items_json}\n")
    else:
        print(f"{config['output_key']}={items_json}")

if __name__ == "__main__":
    main()
