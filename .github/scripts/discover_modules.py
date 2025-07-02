#!/usr/bin/env python
# -*-coding:utf-8 -*-
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
            print(f"Warning: Could not get changed files: {e}")
            print("Running all modules as fallback")
            return None
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
        inputs_file = module_dir / "inputs.json"
        
        if wdl_file.exists() and inputs_file.exists():
            modules.append(module_name)
            print(f"✓ Found valid module: {module_name}")
        else:
            print(f"✗ Skipping {module_name} - missing required files")
            if not wdl_file.exists():
                print(f"    Missing: {wdl_file}")
            if not inputs_file.exists():
                print(f"    Missing: {inputs_file}")
    
    return sorted(modules)

def filter_modules_by_changes(modules, changed_files):
    """Filter modules to only those with changes"""
    if changed_files is None:
        print("No file filtering applied (workflow_dispatch or git error)")
        return modules
        
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

def main():
    print("=== WILDS Module Discovery ===")
    
    # Find all valid modules
    all_modules = find_valid_modules()
    print(f"\nFound {len(all_modules)} valid modules")
    
    # Filter by changes if this is a PR
    changed_files = get_changed_files()
    final_modules = filter_modules_by_changes(all_modules, changed_files)
    
    print(f"\nFinal selection: {len(final_modules)} modules")
    for module in final_modules:
        print(f"  • {module}")
    
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
