#!/usr/bin/env python3
"""Post-process generated documentation HTML files"""

from pathlib import Path
import re
import sys


def fix_badge_paragraphs(html_content: str) -> str:
    """
    Fix badge rendering by removing <p> tags around badge groups.
    The markdown converter wraps badges in <p> tags causing line breaks.
    We want badges to be inline elements with only the <br> controlling line breaks.
    """
    # Pattern to match the opening <p> tag before badges and remove it
    # This looks for <p> followed by badge anchor tags
    pattern = r'<p>(<a href="[^"]*"[^>]*><img src="https://(?:img\.shields\.io|getwilds\.org/badges)[^"]*"[^>]*></a>)'

    # Remove the opening <p> tag before badges
    fixed_content = re.sub(pattern, r'\1', html_content)

    return fixed_content


def postprocess_html_file(file_path: Path) -> None:
    """Post-process a single HTML file."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        original_content = content

        # Apply fixes
        content = fix_badge_paragraphs(content)

        # Only write if content changed
        if content != original_content:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(content)
            print(f"Post-processed: {file_path}")
        else:
            print(f"No changes needed: {file_path}")

    except Exception as e:
        print(f"Error processing {file_path}: {e}", file=sys.stderr)


def main():
    """Post-process all HTML files in the docs directory."""
    docs_dir = Path("docs")

    if not docs_dir.exists():
        print("Error: docs directory not found", file=sys.stderr)
        sys.exit(1)

    html_files = list(docs_dir.glob("**/*.html"))

    if not html_files:
        print("Warning: No HTML files found in docs directory")
        return

    print(f"Found {len(html_files)} HTML files to process")

    for html_file in html_files:
        postprocess_html_file(html_file)

    print("Post-processing complete!")


if __name__ == "__main__":
    main()
