---
name: create-module
description: Create a new WILDS WDL module for a bioinformatics tool following project conventions
argument-hint: <tool-name>
allowed-tools: Bash, Edit, Read, Write, Glob, Grep, WebSearch, WebFetch
---

# Create a New WILDS WDL Module

Create a new module called `ww-$ARGUMENTS` in `modules/ww-$ARGUMENTS/`.

## Steps

### 1. Research the Tool

Before writing any code, research `$ARGUMENTS` to understand:
- What it does and its primary use cases in bioinformatics
- Its command-line interface and key subcommands
- What Docker image exists at `getwilds/$ARGUMENTS` (check Docker Hub or the [wilds-docker-library](https://github.com/getwilds/wilds-docker-library)) and which version tag to pin
- What inputs/outputs are typical

If no `getwilds/` Docker image exists, flag this to the user — one will need to be created first.

If the task would benefit greatly from having an additional tool in the Docker image (e.g., `samtools` for BAM conversion after alignment), stop and let the user know. They can add it to the image in the [wilds-docker-library](https://github.com/getwilds/wilds-docker-library) repo before proceeding.

### 2. Study the Template

Read the canonical template files to match their exact format:
- `modules/ww-template/ww-template.wdl`
- `modules/ww-template/testrun.wdl`
- `modules/ww-template/README.md`

Also review `modules/ww-testdata/ww-testdata.wdl` to see what test data is available (FASTQ, BAM, CRAM, VCF, reference genome, etc.) and pick appropriate test data for this tool.

### 3. Create the Module Files

Create three files in `modules/ww-$ARGUMENTS/`:

#### `ww-$ARGUMENTS.wdl`
- Start with `version 1.0`
- Add a header comment block describing the module
- Define task(s) with the tool's core functionality
- Every task MUST have:
  - `meta` block with: author, email, description, url, outputs
    - Use `author: "WILDS Team"` and `email: "wilds@fredhutch.org"` as placeholders — the actual contributor will replace these in review
  - `parameter_meta` block describing every input
  - `input` block with typed parameters; include `cpu_cores` (Int, default sensible), `memory_gb` (Int, default sensible), and `String docker_image = "getwilds/$ARGUMENTS:<version>"` (pinned, never `latest`)
  - `command <<<` block starting with `set -eo pipefail`, using `~{var}` interpolation
  - `output` block with typed outputs
  - `runtime` block with `docker: docker_image`, `cpu: cpu_cores`, `memory: "~{memory_gb} GB"`

#### `testrun.wdl`
- `version 1.0`
- Import the module using a **relative file path**: `import "ww-$ARGUMENTS.wdl" as ww_$ARGUMENTS_UNDERSCORE`
- Import testdata using a **relative file path**: `import "../ww-testdata/ww-testdata.wdl" as ww_testdata`
- **IMPORTANT**: Use relative imports so the module can be linted and tested locally before being committed. The CI/CD pipeline or the user will update these to GitHub raw URLs when ready.
- Workflow name: `$ARGUMENTS_example` (with hyphens converted to underscores)
- Must run with zero configuration (no input files needed)
- Use appropriate `ww_testdata` tasks to provision test data
- Scatter over samples when applicable
- Define a struct if grouping inputs makes sense

#### `README.md`
- Follow the exact format of `modules/ww-template/README.md`
- Include: badges, overview, module structure, available tasks (with inputs/outputs), usage examples, testing instructions, Docker container info, citations, parameters/resources

### 4. Lint the Module

Run linting and fix any issues. The Makefile requires three tools to be in PATH:
- **sprocket**: `/Users/tfirman/.cargo/bin/sprocket`
- **uv** (used to run miniwdl): `/opt/homebrew/bin/uv`
- **java** (used for WOMtool): `/usr/bin/java`

Always include all three in PATH when running make targets:

```bash
PATH="/Users/tfirman/.cargo/bin:/opt/homebrew/bin:/opt/miniconda3/bin:$PATH" make lint NAME=ww-$ARGUMENTS
```

Fix any lint errors and re-run until clean.

### 5. Test Run (if possible)

Attempt a test run:
```bash
PATH="/Users/tfirman/.cargo/bin:/opt/homebrew/bin:/opt/miniconda3/bin:$PATH" make run_sprocket NAME=ww-$ARGUMENTS
```

If the test run fails due to missing test data or Docker images, note what's needed for the user.

### 6. Important: Do NOT Commit

Do NOT create any git commits or push to GitHub. The user will handle all git operations themselves.

### 7. Summary

Report to the user:
- What files were created
- What Docker image and version was used
- What test data tasks are used
- Any issues that need follow-up (missing Docker image, missing test data in ww-testdata, lint warnings, etc.)
