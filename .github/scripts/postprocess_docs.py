#!/usr/bin/env python3
"""Post-process generated documentation HTML files"""

from pathlib import Path
import re
import sys


def set_default_sidebar_tab(html_content: str) -> str:
    """
    Set the default sidebar tab to 'Full Directory' instead of 'Workflows'.
    Changes showWorkflows from true to false.
    """
    # Find and replace the showWorkflows initialization
    fixed_content = html_content.replace(
        'showWorkflows: $persist(true).using(sessionStorage)',
        'showWorkflows: $persist(false).using(sessionStorage)'
    )

    return fixed_content


def update_page_title(html_content: str, file_path: Path) -> str:
    """
    Update the page title for the homepage from 'Home' to 'WILDS WDL Library'.
    Also updates the homepage header text.
    Only applies to index.html files.
    """
    # Only update title for index.html files
    if file_path.name == 'index.html':
        # Check if this is the root index.html (not in subdirectories)
        if file_path.parent.name == 'docs':
            fixed_content = html_content.replace(
                '<title>Home</title>',
                '<title>WILDS WDL Library</title>'
            )
            # Also update the homepage header
            fixed_content = fixed_content.replace(
                '<h5 class="main__homepage-header">Home</h5>',
                '<h5 class="main__homepage-header">WILDS WDL Library</h5>'
            )
            return fixed_content

    return html_content


def set_accordions_closed(html_content: str) -> str:
    """
    Set all sidebar accordion items to be closed by default, except top-level folders.
    Keeps 'modules' and 'pipelines' visible but collapsed.

    showChildrenCache controls whether folder contents are visible (accordion state).
        - Should be FALSE for all folders (including top-level) so they start collapsed
    showSelfCache controls whether child nodes themselves are shown.
        - Should be TRUE for top-level folders so they're visible
        - Should be FALSE for nested items so they're hidden until parent is expanded
    """
    # Top-level folders that should be visible
    top_level_keys = ["'modules'", "'pipelines'"]

    # Pattern to match showChildrenCache - set ALL to false (all folders start collapsed)
    children_pattern = r"showChildrenCache: \$persist\(\{[^}]+\}\)\.using\(sessionStorage\)"

    def set_all_false(match):
        return match.group(0).replace(': true', ': false')

    fixed_content = re.sub(children_pattern, set_all_false, html_content)

    # Pattern to match showSelfCache - set top-level to true, rest to false
    self_pattern = r"showSelfCache: \$persist\(\{[^}]+\}\)\.using\(sessionStorage\)"

    def set_self_cache(match):
        matched_text = match.group(0)
        # First, set ALL values to false
        result = matched_text.replace(': true', ': false')
        # Then, set top-level folders back to true (they should be visible)
        for key in top_level_keys:
            result = result.replace(f"{key}: false", f"{key}: true")
        return result

    fixed_content = re.sub(self_pattern, set_self_cache, fixed_content)

    return fixed_content


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
        content = update_page_title(content, file_path)
        content = set_default_sidebar_tab(content)
        content = set_accordions_closed(content)
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
