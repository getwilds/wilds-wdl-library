---
name: pr-description
description: Draft a pull request description using the project's PR template
argument-hint: [optional branch or description context]
allowed-tools: Bash, Read, Glob, Grep
---

# Draft a Pull Request Description

Generate a PR description following the project's PR template at `.github/PULL_REQUEST_TEMPLATE.md`.

## Steps

### 1. Gather Context

- Read `.github/PULL_REQUEST_TEMPLATE.md` to get the template format
- Run `git log --oneline main..HEAD` to see all commits on the current branch
- Run `git diff --stat main..HEAD` to understand what files changed
- Optionally read changed files to understand the nature of the changes

### 2. Draft the PR Description

Fill in the PR template sections based on the changes:

- **Type of Change**: Infer from the changes (new module, new pipeline, bug fix, etc.)
- **Description**: Summarize what was added/changed and why, with links to relevant tools/docs where appropriate
- **Related Issue**: Leave as a placeholder unless context is provided
- **Testing**: Describe how the changes were tested based on what you can infer from the conversation or branch history (e.g., sprocket lint, testrun results)
- **Documentation**: Check the boxes based on whether READMEs and parameter descriptions exist in the changed files
- **Additional Context**: Note anything reviewers should know (e.g., Docker image changes, new dependencies)

### 3. Output

Present the draft PR description inside a single markdown code block (```markdown ... ```) so the user can easily copy the raw markdown directly into GitHub. Do NOT render it as formatted text outside a code block. Do NOT create any git commits or push to GitHub.
