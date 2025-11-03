from pathlib import Path
import sys
import re
import os

def is_interactive() -> bool:
    return hasattr(sys, "ps1")

def get_root() -> Path:
    if is_interactive():
        return Path.cwd().resolve()
    else:
        return Path(__file__).parent.parent.parent.resolve()

def fetch_file(path: Path) -> list[str]:
    try:
        with open(path, "r", encoding="utf-8") as f:
            return f.readlines()
    except FileNotFoundError:
        return []


def convert_github_imports_to_relative(wdl_content: list[str], wdl_path: Path, root: Path) -> list[str]:
    """Convert GitHub URL imports to relative paths."""
    # Pattern to match GitHub imports - capture everything after the branch name
    github_pattern = r'import\s+"https://raw\.githubusercontent\.com/getwilds/wilds-wdl-library/refs/heads/[^/]+/(.+?)"\s+as\s+'

    converted_lines = []
    for line in wdl_content:
        match = re.search(github_pattern, line)
        if match:
            # Extract the path portion (e.g., "modules/ww-testdata/ww-testdata.wdl")
            github_path = match.group(1)

            # Calculate relative path from current WDL file to the imported file
            current_dir = wdl_path.parent
            target_path = root / github_path

            # Use os.path.relpath to compute the relative path correctly
            relative_path = os.path.relpath(target_path, current_dir)

            # Replace the GitHub URL with relative path
            new_line = re.sub(
                r'import\s+"https://raw\.githubusercontent\.com/getwilds/wilds-wdl-library/refs/heads/[^/]+/.+?"\s+as\s+',
                f'import "{relative_path}" as ',
                line
            )
            converted_lines.append(new_line)
        else:
            converted_lines.append(line)

    return converted_lines

def add_docs_callout_to_readme(readme_path: Path) -> None:
    """Add documentation-specific callout to the main README."""
    readme = fetch_file(readme_path)

    # Find where to insert the callout (after the badges section, before ## Overview)
    overview_index = None
    for i, line in enumerate(readme):
        if line.strip().startswith("## Overview"):
            overview_index = i
            break

    if overview_index is None:
        return  # Can't find insertion point, skip

    # Create the callout box
    callout = [
        "\n",
        "---\n",
        "\n",
        "> **📚 Viewing the Documentation Site?**\n",
        ">\n",
        "> Welcome to the WILDS WDL Library technical documentation! This site provides comprehensive documentation for all available WDL modules, vignettes, and workflows.\n",
        ">\n",
        "> **Getting Started:**\n",
        "> - Browse available components using the **sidebar navigation**\n",
        "> - Each module page includes task descriptions, inputs/outputs, and usage examples\n",
        "> - **Vignettes** demonstrate how to combine modules into common analysis patterns\n",
        "> - **Workflows** provide complete, production-ready analysis pipelines\n",
        ">\n",
        "> **Using Components in Your Workflows:**\n",
        "> - Import modules directly via GitHub URLs (see examples below)\n",
        "> - All components are tested with Cromwell, miniWDL, and Sprocket\n",
        "> - For contributing or development setup, see the sections below\n",
        "\n",
        "---\n",
        "\n",
    ]

    # Insert the callout before ## Overview
    readme = readme[:overview_index] + callout + readme[overview_index:]

    # Write the modified README
    with open(readme_path, "w", encoding="utf-8") as f:
        f.writelines(readme)


def make_preambles() -> None:
    root = get_root()

    # First, modify the main README to add docs-specific callout
    main_readme = root / "README.md"
    if main_readme.exists():
        add_docs_callout_to_readme(main_readme)

    # Then process all WDL files
    paths = list(root.glob("**/*.wdl"))

    for path in paths:
        # Skip testrun.wdl files
        if path.name == "testrun.wdl":
            continue

        path_readme = list(path.parent.glob("README.md"))
        readme = fetch_file(path_readme[0])
        readme_prefixed = [f"## {w}" for w in readme]

        wdl = fetch_file(path)

        # Convert GitHub imports to relative paths
        wdl_converted = convert_github_imports_to_relative(wdl, path, root)

        readme_wdl = readme_prefixed + ["\n", "\n"] + wdl_converted

        with open(path, "w", encoding="utf-8") as f:
            f.writelines(readme_wdl)


if __name__ == "__main__":
    make_preambles()
