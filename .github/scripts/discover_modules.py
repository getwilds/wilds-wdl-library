#!/usr/bin/env -S uv run --script
#
# /// script
# requires-python = "==3.13"
# ///
'''
@File    :   discover_modules.py
@Time    :   2025/07/01 22:47:23
@Author  :   Taylor Firman
@Version :   v0.1
@Contact :   tfirman@fredhutch.org
@Desc    :   Identifies WILDS modules needing a test run based on changes in a PR.
'''

import json
import os
import subprocess
import sys
from pathlib import Path

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

def find_valid_modules():
    """Find all directories with valid module structure"""
    modules = []
    modules_dir = Path('modules')
    
    if not modules_dir.exists():
        print("No modules directory found")
        return modules
        
    for module_dir in modules_dir.iterdir():
        if not module_dir.is_dir():
            continue
            
        module_name = module_dir.name
        wdl_file = module_dir / f"{module_name}.wdl"
        
        if wdl_file.exists():
            modules.append(module_name)
            print(f"Found valid module: {module_name}")
        else:
            print(f"Skipping {module_name} - missing required files")
            if not wdl_file.exists():
                print(f"    Missing: {wdl_file}")
    
    return sorted(modules)

def filter_modules_by_changes(modules, changed_files):
    """Filter modules to only those with changes"""
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
    for module in modules:
        module_files = [f for f in changed_files if f.startswith(f'modules/{module}/')]
        if module_files:
            filtered.append(module)
            print(f"Including {module} ({len(module_files)} changed files)")
        else:
            print(f"Skipping {module} (no changes)")
            
    return filtered

def handle_workflow_dispatch_input(modules, dispatch_module=None):
    """Handle manual workflow dispatch input for specific module"""
    if dispatch_module is None:
        dispatch_module = ''
    dispatch_module = dispatch_module.strip()
    
    if not dispatch_module:
        print("No specific module requested for workflow_dispatch, running all modules")
        return modules
    
    # Check if the requested module exists and is valid
    if dispatch_module in modules:
        print(f"Running specific module from workflow_dispatch: {dispatch_module}")
        return [dispatch_module]
    else:
        print(f"ERROR: Requested module '{dispatch_module}' not found or invalid")
        print(f"Available modules: {', '.join(modules)}")
        # Return empty list to fail gracefully rather than running wrong modules
        return []

def main():
    print("=== WILDS Module Discovery ===")
    
    # Parse command line arguments
    dispatch_module = None
    if len(sys.argv) > 1:
        dispatch_module = sys.argv[1]
        print(f"Workflow dispatch module argument: '{dispatch_module}'")
    
    # Find all valid modules
    all_modules = find_valid_modules()
    print(f"\nFound {len(all_modules)} valid modules")
    
    # Check event type
    event_name = os.environ.get('GITHUB_EVENT_NAME', '')
    print(f"GitHub event: {event_name}")
    
    if event_name == 'workflow_dispatch':
        # Handle manual workflow dispatch
        final_modules = handle_workflow_dispatch_input(all_modules, dispatch_module)
    else:
        # Handle PR - filter by changes
        changed_files = get_changed_files()
        final_modules = filter_modules_by_changes(all_modules, changed_files)
    
    print(f"\nFinal selection: {len(final_modules)} modules")
    for module in final_modules:
        print(f"  - {module}")
    
    # Output for GitHub Actions
    modules_json = json.dumps(final_modules)
    
    # Write to GITHUB_OUTPUT
    github_output = os.environ.get('GITHUB_OUTPUT')
    if github_output:
        with open(github_output, 'a') as f:
            f.write(f"modules={modules_json}\n")
    else:
        print(f"modules={modules_json}")

if __name__ == "__main__":
    main()
