# AGENTS.md - WILDS WDL Library

## Project Overview
Centralized bioinformatics WDL library with two tiers: **modules** (tool-specific task collections) and **pipelines** (workflows combining modules). All WDL uses version 1.0.

## Repo Structure
```
modules/ww-toolname/
  ww-toolname.wdl          # Task definitions
  testrun.wdl              # Zero-config test workflow (required)
  README.md
pipelines/ww-pipeline-name/
  ww-pipeline-name.wdl     # Workflow importing modules
  testrun.wdl
  inputs.json              # Example inputs with placeholder paths
  README.md
```

## Module Conventions
- Use `ww-template` as the canonical example
- Every task needs `meta` (author, email, description, url, outputs) and `parameter_meta` blocks
- Command blocks: start with `set -eo pipefail`, use `<<<...>>>` heredoc syntax
- Docker images: always from `getwilds/` with exact version pins (never `latest`)
- Resource params: `cpu_cores` and `memory_gb` with sensible defaults
- Imports use GitHub raw URLs: `https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-*/ww-*.wdl`
- Modules are self-contained: `ww-<toolname>.wdl` files never import other modules. Cross-module composition happens only in `testrun.wdl` (for test fixtures) and in pipelines. If a task needs another tool's output, add it to the existing module rather than introducing a cross-module import in the source WDL.

## Test Workflows (testrun.wdl)
- Must run with zero configuration (no input files needed)
- Import the module under test + `ww-testdata` for auto-provisioned test data
- Workflow name: `{toolname}_example` (e.g., `star_example`)
- Scatter over samples when applicable
- Some modules also provide a `testrun_hpc.wdl` variant for tools that require HPC-specific resources (GPUs, large memory)

## Pipeline Conventions
- Combine existing modules; prefer existing modules over new task definitions
- Document complexity level: Basic (2-3 modules), Intermediate (4-6), Advanced (10+)
- Include `inputs.json` with placeholder paths as user template
- Use WDL structs for grouping complex inputs

## Linting & Testing
```bash
make lint NAME=ww-toolname          # All linting (Sprocket + miniwdl + WOMtool + Cirro)
make lint_sprocket NAME=ww-toolname # Sprocket only
make run NAME=ww-toolname           # Run testrun.wdl on all engines
make run_sprocket NAME=ww-toolname  # Run testrun.wdl with Sprocket
make run_miniwdl NAME=ww-toolname   # Run testrun.wdl with miniwdl
make run_cromwell NAME=ww-toolname  # Run testrun.wdl with Cromwell
```
Sprocket exceptions (from sprocket.toml): TodoComment, ContainerUri, UnusedInput

## Common Development Pitfalls
- Import URLs must point to the correct branch during development; switch to `main` before merging
- Docker image version conflicts (especially with complex tools like ColabFold) - always pin versions
- Test data changes go in the `ww-testdata` module, not individual modules
- All modules/pipelines must pass CI on three executors: Cromwell, miniWDL, and Sprocket

## Dockstore Registration
- All modules and pipelines must be registered in `.dockstore.yml`
- Modules go under `tools:`, pipelines under `workflows:`
- Include author info (name, email, role, affiliation, orcid), description, topic, and tags
- Keep entries alphabetically ordered within their section
- Author details should match `CITATION.cff` where applicable

## PR Process
- Request reviews from @emjbishop or @tefirman
- Fill out PR template (type, description, testing, documentation checklist)
- CI runs linting + multi-executor test runs automatically

## Agent Skills

Reusable task-specific instructions live in [.agents/skills/](.agents/skills/). Each `SKILL.md` is a self-contained recipe an agent can follow:

- `create-module` — scaffold a new `ww-toolname` module
- `create-pipeline` — assemble existing modules into a new pipeline
- `add-testdata` — add a test-data download task to `ww-testdata`
- `run-tests` — run `testrun.wdl` via sprocket or miniwdl
- `lint-module` — run linting and fix issues
- `pr-description` — draft a PR description from the current branch
