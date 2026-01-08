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


def add_github_link(html_content: str) -> str:
    """
    Add a GitHub icon link just to the left of the theme toggle button.
    Links to the wilds-wdl-library repository.
    """
    # Skip if already added
    if 'title="View on GitHub"' in html_content:
        return html_content

    # GitHub icon styled to match the theme toggle button
    github_icon = '''<a href="https://github.com/getwilds/wilds-wdl-library" target="_blank" rel="noopener noreferrer" class="border border-slate-700 rounded-md h-8 flex items-center justify-center text-slate-300 w-8 cursor-pointer hover:border-slate-500" style="margin-right: 8px;" title="View on GitHub"><svg fill="currentColor" viewBox="0 0 24 24" style="width: 1.25rem; height: 1.25rem;"><path d="M12 0c-6.626 0-12 5.373-12 12 0 5.302 3.438 9.8 8.207 11.387.599.111.793-.261.793-.577v-2.234c-3.338.726-4.033-1.416-4.033-1.416-.546-1.387-1.333-1.756-1.333-1.756-1.089-.745.083-.729.083-.729 1.205.084 1.839 1.237 1.839 1.237 1.07 1.834 2.807 1.304 3.492.997.107-.775.418-1.305.762-1.604-2.665-.305-5.467-1.334-5.467-5.931 0-1.311.469-2.381 1.236-3.221-.124-.303-.535-1.524.117-3.176 0 0 1.008-.322 3.301 1.23.957-.266 1.983-.399 3.003-.404 1.02.005 2.047.138 3.006.404 2.291-1.552 3.297-1.23 3.297-1.23.653 1.653.242 2.874.118 3.176.77.84 1.235 1.911 1.235 3.221 0 4.609-2.807 5.624-5.479 5.921.43.372.823 1.102.823 2.222v3.293c0 .319.192.694.801.576 4.765-1.589 8.199-6.086 8.199-11.386 0-6.627-5.373-12-12-12z"/></svg></a>'''

    # Wrap the theme toggle button and add GitHub icon before it, in a flex container
    # This keeps both buttons together as a single flex item
    html_content = re.sub(
        r'''(<button x-on:click="\s*document\.documentElement\.classList\.toggle\('light'\)\s*localStorage\.setItem\('theme', document\.documentElement\.classList\.contains\('light'\) \? 'light' : 'dark'\)\s*" class="[^"]*">☀︎</button>)''',
        rf'<div style="display: flex; align-items: center;">{github_icon}\1</div>',
        html_content
    )

    return html_content


def shorten_sprocket_paths(html_content: str) -> str:
    """
    Shorten full WDL paths in sprocket run commands to just the filename.
    Only affects paths within the 'RUN WITH' section.
    e.g., 'modules/ww-annotsv/ww-annotsv.wdl' becomes 'ww-annotsv.wdl'
    """
    # Match paths within the run-with-content section
    # The path is inside a <span> tag after "sprocket run"
    # Unix paths: modules/name/file.wdl or pipelines/name/file.wdl
    html_content = re.sub(
        r'(main__run-with-content.*?)(?:modules|pipelines)/[^/]+/([^/<]+\.wdl)',
        r'\1\2',
        html_content
    )
    # Windows paths: modules\name\file.wdl or pipelines\name\file.wdl
    html_content = re.sub(
        r'(main__run-with-content.*?)(?:modules|pipelines)\\[^\\]+\\([^\\<]+\.wdl)',
        r'\1\2',
        html_content
    )
    return html_content


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
        content = shorten_sprocket_paths(content)
        content = add_github_link(content)

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
