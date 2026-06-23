---
name: lint-module
description: Run WDL linting on a WILDS module or pipeline and fix any issues found
argument-hint: <ww-name>
allowed-tools: Bash(make *), Read, Edit, Glob, Grep
---

# Lint a WILDS WDL Module or Pipeline

Run all available linting checks on `$ARGUMENTS` and fix any issues.

## Steps

### 1. Determine Module or Pipeline

Check if `$ARGUMENTS` exists in `modules/` or `pipelines/`:
- Module: `modules/$ARGUMENTS/$ARGUMENTS.wdl`
- Pipeline: `pipelines/$ARGUMENTS/$ARGUMENTS.wdl`

### 2. Run Linting

```bash
make lint NAME=$ARGUMENTS
```

This runs sprocket, miniwdl, and WOMtool linters.

If make targets aren't available, run sprocket directly:
```bash
make lint_sprocket NAME=$ARGUMENTS
```

### 3. Fix Issues

For each lint error:
- Read the relevant WDL file
- Fix the issue following project conventions (see CLAUDE.md)
- Re-run linting to confirm the fix

**Known Sprocket exceptions** (these are OK and won't cause failures):
- TodoComment
- ContainerUri
- TrailingComma
- CommentWhitespace
- UnusedInput

### 4. Report Results

Summarize what was found and fixed. If there are unfixable issues, explain why and suggest next steps.
