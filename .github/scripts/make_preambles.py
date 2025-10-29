from pathlib import Path
import sys
import re

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
    # Pattern to match GitHub imports
    github_pattern = r'import\s+"https://raw\.githubusercontent\.com/getwilds/wilds-wdl-library/[^"]+/([^"]+)"\s+as\s+'

    converted_lines = []
    for line in wdl_content:
        match = re.search(github_pattern, line)
        if match:
            # Extract the path portion (e.g., "modules/ww-testdata/ww-testdata.wdl")
            github_path = match.group(1)

            # Calculate relative path from current WDL file to the imported file
            current_dir = wdl_path.parent
            target_path = root / github_path

            try:
                relative_path = target_path.relative_to(current_dir)
            except ValueError:
                # If relative_to fails, use a fallback approach
                relative_path = Path("..") / github_path

            # Replace the GitHub URL with relative path
            new_line = re.sub(
                r'import\s+"https://raw\.githubusercontent\.com/getwilds/wilds-wdl-library/[^"]+/[^"]+"\s+as\s+',
                f'import "{relative_path}" as ',
                line
            )
            converted_lines.append(new_line)
        else:
            converted_lines.append(line)

    return converted_lines

def make_preambles() -> None:
    root = get_root()
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
