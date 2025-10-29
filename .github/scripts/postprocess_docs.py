#!/usr/bin/env python3
"""Post-process generated documentation HTML files"""

from pathlib import Path
import re
import sys


def fix_badge_paragraphs(html_content: str) -> str:
    """
    Fix badge rendering by wrapping all badges in a single div with inline display.
    The markdown converter creates line breaks between badges, we want them inline.
    """
    # Pattern to match the entire badge section from first badge to closing </p>
    # This captures: optional <p>, all badge links (with <br>), and closing </p>
    pattern = r'(?:<p>)?(<a href="[^"]*"[^>]*><img src="https://(?:img\.shields\.io|getwilds\.org/badges)[^"]*"[^>]*></a>(?:\s*\n\s*<a href="[^"]*"[^>]*><img src="https://(?:img\.shields\.io|getwilds\.org/badges|github\.com/[^/]+/[^/]+/actions)[^"]*"[^>]*></a>|\s*<br>\s*)*)</p>'

    def replace_badges(match):
        badge_content = match.group(1)
        # Remove newlines between badges (but keep <br> tags)
        badge_content = re.sub(r'>\s*\n\s*<a', '> <a', badge_content)

        # Add inline display style to each anchor tag to force inline rendering
        badge_content = re.sub(
            r'<a href',
            r'<a style="display: inline-block; margin-right: 4px;" href',
            badge_content
        )

        # Wrap in a div with inline styling to keep badges inline
        return f'<div style="line-height: 2.5; margin-bottom: 1em;">{badge_content}</div>'

    # Replace the badge paragraph with a styled div
    fixed_content = re.sub(pattern, replace_badges, html_content, flags=re.MULTILINE | re.DOTALL)

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
