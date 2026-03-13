#!/usr/bin/env python3
"""
Validate .cirro configurations in WILDS pipelines.

Checks that .cirro directories contain the required files,
JSON files are valid, and preprocess.py has no syntax errors.
"""

import json
import sys
from pathlib import Path

REQUIRED_FILES = [
    "preprocess.py",
    "process-form.json",
    "process-input.json",
    "process-output.json",
    "process-compute.config",
]

def validate_json_file(filepath):
    """Validate that a file contains valid JSON. Returns list of error strings."""
    errors = []
    try:
        with open(filepath) as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        errors.append(f"  Invalid JSON in {filepath.name}: {e}")
        return errors, None
    return errors, data


def validate_form(filepath):
    """Validate process-form.json has expected structure."""
    errors, data = validate_json_file(filepath)
    if data is None:
        return errors

    if not isinstance(data, dict):
        errors.append(f"  {filepath.name}: expected a JSON object at top level")
        return errors

    if "form" not in data:
        errors.append(f"  {filepath.name}: missing top-level 'form' key")
        return errors

    form = data["form"]
    if not isinstance(form, dict):
        errors.append(f"  {filepath.name}: 'form' should be an object")
        return errors

    if "properties" not in form:
        errors.append(f"  {filepath.name}: 'form' missing 'properties' key")

    if "required" in form and not isinstance(form["required"], list):
        errors.append(f"  {filepath.name}: 'required' should be a list")

    return errors


def validate_input(filepath):
    """Validate process-input.json has JSON path mappings."""
    errors, data = validate_json_file(filepath)
    if data is None:
        return errors

    if not isinstance(data, dict):
        errors.append(f"  {filepath.name}: expected a JSON object")
        return errors

    for key, value in data.items():
        if not isinstance(value, str):
            errors.append(f"  {filepath.name}: value for '{key}' should be a string, got {type(value).__name__}")
        elif not value.startswith("$."):
            errors.append(f"  {filepath.name}: value for '{key}' should be a JSON path (start with '$.')")

    return errors


def validate_output(filepath):
    """Validate process-output.json is valid JSON."""
    errors, _ = validate_json_file(filepath)
    return errors


def validate_preprocess(filepath):
    """Validate preprocess.py has no syntax errors."""
    errors = []
    try:
        source = filepath.read_text()
        compile(source, str(filepath), "exec")
    except SyntaxError as e:
        errors.append(f"  {filepath.name}: Python syntax error: {e}")
    return errors


def validate_cirro_dir(cirro_dir):
    """Validate a single .cirro directory. Returns list of error strings."""
    errors = []

    # Check required files
    for filename in REQUIRED_FILES:
        if not (cirro_dir / filename).exists():
            errors.append(f"  Missing required file: {filename}")

    # Validate individual files
    form_path = cirro_dir / "process-form.json"
    if form_path.exists():
        errors.extend(validate_form(form_path))

    input_path = cirro_dir / "process-input.json"
    if input_path.exists():
        errors.extend(validate_input(input_path))

    output_path = cirro_dir / "process-output.json"
    if output_path.exists():
        errors.extend(validate_output(output_path))

    preprocess_path = cirro_dir / "preprocess.py"
    if preprocess_path.exists():
        errors.extend(validate_preprocess(preprocess_path))

    return errors


def main():
    pipelines_dir = Path("pipelines")
    if not pipelines_dir.exists():
        print("No pipelines directory found")
        return 0

    found_any = False
    all_errors = {}

    for pipeline_dir in sorted(pipelines_dir.iterdir()):
        if not pipeline_dir.is_dir():
            continue

        cirro_dir = pipeline_dir / ".cirro"
        if not cirro_dir.is_dir():
            print(f"Skipping {pipeline_dir.name} (no .cirro directory)")
            continue

        found_any = True
        print(f"Validating {pipeline_dir.name}/.cirro/ ...")
        errors = validate_cirro_dir(cirro_dir)

        if errors:
            all_errors[pipeline_dir.name] = errors
            print(f"  FAIL ({len(errors)} issue(s))")
        else:
            print(f"  OK")

    if not found_any:
        print("No .cirro directories found in any pipeline")
        return 0

    if all_errors:
        print(f"\n{'='*50}")
        print(f"Cirro validation failed for {len(all_errors)} pipeline(s):\n")
        for pipeline, errors in all_errors.items():
            print(f"{pipeline}:")
            for error in errors:
                print(error)
            print()
        return 1

    print(f"\nAll Cirro configurations valid!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
