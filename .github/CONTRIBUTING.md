# Contributing to the WILDS WDL Library

Thank you for your interest in contributing to the WILDS WDL Library! This document provides guidelines for contributing modules, vignettes, workflows, and improvements to our centralized collection of bioinformatics WDL infrastructure.

## Table of Contents

- [Getting Started](#getting-started)
- [Repository Structure](#repository-structure)
- [Types of Contributions](#types-of-contributions)
- [Module Development Guidelines](#module-development-guidelines)
- [Vignette Development Guidelines](#vignette-development-guidelines)
- [Workflow Development Guidelines](#workflow-development-guidelines)
- [Testing Requirements](#testing-requirements)
- [Documentation Standards](#documentation-standards)
- [Documentation Website](#documentation-website)
- [Pull Request Process](#pull-request-process)
- [Code of Conduct](#code-of-conduct)

## Getting Started

Before contributing code changes, please:

1. **[Fork the repository](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo)** to your GitHub account

2. **Set up your development environment** with the required tools:
   - For local testing:
       - [sprocket](https://sprocket.bio/installation.html) (recommended)
       - [miniWDL](https://miniwdl.readthedocs.io/en/latest/getting_started.html#install-miniwdl)
       - [uv](https://docs.astral.sh/uv/getting-started/installation/) for automated local testing
   - [Docker Desktop](https://www.docker.com/get-started/) for container execution

3. **Make code changes** and push them to your fork

4. **Submit a pull request (PR)** to merge your contributions into the `main` branch of the original repo
    *   The title of your PR should briefly describe the change.
    *   If your contribution resolves an issue, the body of your PR should contain `Fixes #issue-number`

## Repository Structure

The WILDS WDL Library follows a three-tier architecture:

- **Modules**: Collection of tasks that use a given tool
- **Vignettes**: Small example pipelines that import module tasks
- **Workflows**: Production-ready pipelines

```
wilds-wdl-library/
├── modules/
│   └── ww-toolname/
│       ├── ww-toolname.wdl
│       └── README.md
├── vignettes/
│   └── ww-example-name/
│       ├── ww-example-name.wdl
│       ├── inputs.json
│       └── README.md
├── workflows/
│   └── ww-pipeline-name/
│       ├── ww-pipeline-name.wdl
│       ├── inputs.json
│       └── README.md
└── .github/
    └── workflows/     # CI/CD automation
```

## Types of Contributions

### 1. Bug Reports and Issues

- Use the [GitHub Issues](https://github.com/getwilds/wilds-wdl-library/issues) page
- Provide detailed information about the problem
- Include error messages, info about input files, and steps to reproduce
- Tag issues appropriately (bug, enhancement, question, etc.)

### 2. Documentation Improvements

- Fix typos, improve clarity, or add missing information
- Enhance README files with better examples

### 3. Module Contributions

- Focus on one high-utility bioinformatics tool
- Follow standardized module structure
- Include comprehensive testing and validation

### 4. Vignette Contributions

- Combine 2-3 existing modules
- Demonstrate common analysis patterns
- Provide educational value for the community

### 5. Workflow Contributions

- Create publication-ready analysis pipelines
- Integrate multiple modules
- Include comprehensive validation with realistic datasets

## Module Development Guidelines

**See our [ww-template module](https://github.com/getwilds/wilds-wdl-library/tree/main/modules/ww-template) as an example**

**The module folder must contain:**

1. **`ww-toolname.wdl`** - Main WDL file containing task definitions for the tool
2. **`testrun.wdl`** - Test workflow demonstrating module functionality (must be named `testrun.wdl`)
3. **`README.md`** - Comprehensive documentation


**Your main WDL file (`ww-toolname.wdl`) must include:**

- **Version declaration**: Use WDL version 1.0
- **Task definitions**: Individual tasks with proper resource requirements
- **Metadata documentation**: Describe properties of tasks (e.g. inputs, outputs) using `meta` and `parameter_meta` blocks

**Your test workflow file (`testrun.wdl`) must include:**

- **Version declaration**: Use WDL version 1.0
- **Module imports**: Import the module being tested and the `ww-testdata` module using GitHub URLs
- **Sample struct definition**: Define a struct for organizing sample inputs if needed
- **Test workflow**: A `toolname_example` workflow that demonstrates all tasks (must follow the naming convention `{module}_example` where `{module}` is the tool name, e.g., `star_example` for `ww-star`)
- **Auto-downloading of test data**: Use the `ww-testdata` module to automatically provision test data
- **Validation task (optional)**: Consider including a validation task to verify output correctness

**Parameter preferences:**

- Use descriptive parameter names
- Include optional parameters with sensible defaults
- Support both single samples and batch processing where applicable

**Docker image preferences:**

- Use images from the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library) when available
- If creating new images, follow WILDS container standards and consider contributing to the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library).
- Specify exact image versions (avoid `latest` tags)
- Document image dependencies in the README


## Vignette Development Guidelines

**Vignettes should:**

- Combine 2-3 existing modules from the library
- Demonstrate realistic analysis workflows
- Serve as educational templates
- Use publicly available test data

**No New Tasks Rule**

- Vignettes should **only** combine existing modules - do not create new task definitions. If you need new functionality, contribute it as a module first.

**Vignette inputs.json**

In order to pass the automated GitHub Action (GHA) tests, your `inputs.json` must:

- Use no more than 4 CPUs (for the sprocket executor)
- Use input files that get downloaded as part of the GHA.
    - See the `download-test-data` job within `./github/workflows/vignettes-testrun.yml`
    - Modify `vignettes-testrun.yml` as needed

**Platform-Specific Configurations (Optional)**

Vignettes may include optional platform-specific configuration directories for execution on cloud platforms or workflow management systems:

- **Location**: Place platform configs in a subdirectory within the vignette (e.g., `vignettes/ww-example/.cirro/`)
- **Naming convention**: Use dotfile directory names (`.cirro/`, `.terra/`, etc.) to indicate platform
- **Standalone principle**: Keep all vignette-related files (WDL, inputs, platform configs) in the vignette directory
- **Documentation**: Document platform configurations in the vignette's README with links to platform documentation
- **Examples**:
  - `.cirro/` for [Cirro](https://cirro.bio/) platform ([config documentation](https://docs.cirro.bio/pipelines/adding-pipelines/))
  - `.terra/` for [Terra](https://terra.bio/) workspace configurations
  - Other platform-specific directories as needed

Platform configurations are entirely optional and should not be required to run the vignette with standard WDL executors (Cromwell, miniWDL, Sprocket).

## Workflow Development Guidelines

**Workflows represent the highest tier and must meet production standards:**

- **Comprehensive validation**: Extensive testing with realistic datasets
- **Scalability**: Support for large datasets and batch processing
- **Documentation**: Thorough documentation to aid in publication

## Testing Requirements

### Local Tests

Make sure you have these installed:

- [sprocket](https://sprocket.bio/installation.html) (recommended)
- [miniWDL](https://miniwdl.readthedocs.io/en/latest/getting_started.html#install-miniwdl)
- [uv](https://docs.astral.sh/uv/getting-started/installation/) for automated testing with our Makefile
- [Docker Desktop](https://www.docker.com/get-started/) for container execution

#### Option 1: Manual Testing

Test your WDL manually by navigating to the module directory:

```bash
cd modules/ww-toolname

# Linting with miniwdl (check both main module and test workflow)
miniwdl check ww-toolname.wdl
miniwdl check testrun.wdl

# Linting with sprocket (ignoring things we don't care about)
sprocket lint \
  -e TodoComment \
  -e ContainerUri \
  -e TrailingComma \
  -e CommentWhitespace \
  -e UnusedInput \
  ww-toolname.wdl

sprocket lint \
  -e TodoComment \
  -e ContainerUri \
  -e TrailingComma \
  -e CommentWhitespace \
  -e UnusedInput \
  testrun.wdl

# Test running (use testrun.wdl for execution tests)
sprocket run testrun.wdl --entrypoint toolname_example
miniwdl run testrun.wdl
```

#### Option 2: Automated Testing with Makefile (Recommended)

Use our automated Makefile from the repository root for easier testing:

```bash
# Test a specific module (replace ww-toolname with your module name)
make lint MODULE=ww-toolname          # Run all linting checks
make lint_sprocket MODULE=ww-toolname # Run only sprocket linting
make lint_miniwdl MODULE=ww-toolname  # Run only miniwdl linting
make run_sprocket MODULE=ww-toolname  # Run sprocket with proper entrypoint
make run_miniwdl MODULE=ww-toolname   # Run miniwdl

# Test all modules
make lint    # Lint all modules
make run     # Run all modules with both sprocket and miniwdl
```

The Makefile automatically handles:
- Proper entrypoint naming for sprocket (`{module}_example`)
- Module discovery and validation
- Dependency checking (sprocket, uv, etc.)
- Consistent test execution across all modules

### Test Data

- Use the [`ww-testdata` module](https://github.com/getwilds/wilds-wdl-library/tree/main/modules/ww-testdata) for standardized test datasets
- If you need additional test datasets, modify the `ww-testdata` module also
- Include small, representative test files in your examples

### Automated Tests

All contributions must pass our automated testing pipeline which executes on a PR via GitHub Actions:

- **Multi-executor validation**: Tests with Cromwell, miniWDL, and Sprocket
- **Container verification**: All Docker images must be accessible and functional
- **Syntax validation**: WDL syntax and structure validation
- **Integration testing**: Cross-module compatibility testing

## Documentation Website

The WILDS WDL Library includes an automatically-generated documentation website that provides comprehensive technical documentation for all modules, vignettes, and workflows. Understanding how this documentation works is important for contributors.

### How Documentation is Generated

The documentation website is built using [Sprocket](https://sprocket.bio/) and automatically deployed to GitHub Pages. The documentation is generated from:

- **README files**: Each module, vignette, and workflow directory contains a README.md that becomes the documentation homepage for that component
- **WDL files**: Task descriptions, inputs, outputs, and metadata are automatically extracted from WDL files
- **Main README**: The repository's root README.md serves as the documentation site homepage

### Automatic Deployment

Documentation is automatically built and deployed when changes are merged to the `main` branch:

1. The [build-docs.yml](.github/workflows/build-docs.yml) GitHub Actions workflow triggers on push to `main`
2. The workflow runs the [make_preambles.py](.github/scripts/make_preambles.py) script to prepare WDL files
3. Sprocket generates static HTML documentation
4. The [postprocess_docs.py](.github/scripts/postprocess_docs.py) script applies final formatting
5. Documentation is deployed to GitHub Pages at the repository's documentation URL

**Important**: You don't need to build or commit documentation files - they are generated automatically in CI/CD.

### Previewing Documentation Locally

Before submitting a PR, you can preview how your changes will appear on the documentation website using the provided Makefile targets:

#### Build and Preview Documentation

```bash
# Build documentation locally (mirrors the CI/CD process)
make docs-preview

# Serve the documentation on http://localhost:8000
make docs-serve

# Or do both in one command
make docs
```

The `docs-preview` target will:
- Check for uncommitted changes and warn you (docs are built from your last commit)
- Safely stash any uncommitted work
- Run the same build process as the GitHub Actions workflow
- Generate documentation in the `docs/` directory
- Restore your uncommitted changes when finished
- Clean up all temporary build files

**Note**: The `docs/` directory is gitignored and should never be committed to the repository.

#### What Gets Built

When you run `make docs-preview`, the build process:
1. Prepends each module's README to its WDL file for better documentation context
2. Converts GitHub import URLs to relative paths for local navigation
3. Generates comprehensive HTML documentation for all tasks, workflows, and components
4. Applies custom styling and post-processing

### Documentation Best Practices

When contributing, ensure your documentation is clear and complete:

- **README files**: Write clear, user-focused descriptions of what your module/vignette/workflow does
- **Task metadata**: Use `meta` blocks to document task purpose, authors, and other high-level information
- **Parameter metadata**: Use `parameter_meta` blocks to describe all inputs and outputs
- **Examples**: Include usage examples in README files
- **Preview locally**: Always run `make docs-preview` before submitting a PR to verify how your documentation will appear

### Troubleshooting Documentation Builds

If you encounter issues with local documentation builds:

- Ensure you have the required dependencies installed (`sprocket`, `uv`, `python 3.13`)
- Check that you're running the command from the repository root
- Review error messages - they often indicate issues with WDL syntax or README formatting

For questions about documentation, please contact [wilds@fredhutch.org](mailto:wilds@fredhutch.org).

## Pull Request Process

After meeting the requirements above, submit a PR to merge your forked repo into `main`.

1. **Create descriptive PR title**:
   - Examples: `Add BWA alignment module`, `RNA-seq analysis example vingette`

2. **Fill out PR template**: Provide detailed information about your contribution

3. **Link related issues**: Reference any GitHub issues your PR addresses

4. **Request reviews**: Tag Emma Bishop (@emjbishop) or Taylor Firman (@tefirman)

### Review Criteria

Your PR will be evaluated on:

- **Functionality**: Does it work as intended?
- **Testing**: Are tests comprehensive and passing?
- **Documentation**: Is documentation clear and complete?
- **Standards compliance**: Does it follow WILDS conventions?
- **Code quality**: Is the WDL code well-structured and readable?
- **Uniqueness**: Does it avoid duplicating existing functionality in the library?

## Help for new contributors

New contributors are welcome! If you're new to WDL or bioinformatics workflows:

- Review our [WDL 101 course materials](https://github.com/getwilds/wdl-101)
- Check out existing modules for examples
- Don't hesitate to ask questions in issues or via email. If you have a `uw.edu` or `fredhutch.org` email you can also ask questions in our `fh-data` [slack workspace]([https://hutchdatascience.org/joinslack/)
- Consider starting with documentation contributions

For more questions you can contact the Fred Hutch Data Science Lab at [wilds@fredhutch.org](mailto:wilds@fredhutch.org)


## Code of Conduct

By participating in this project, you agree to abide by our [code of conduct](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CODE_OF_CONDUCT.md):

- **Be respectful**: Treat all community members with respect and kindness
- **Be collaborative**: Work together constructively and help others learn
- **Be inclusive**: Welcome contributors from all backgrounds and experience levels
- **Be patient**: Remember that everyone is learning and growing

### Reporting Issues

If you experience or witness unacceptable behavior, please report it to [wilds@fredhutch.org](mailto:wilds@fredhutch.org).

## License

By contributing to this project, you agree that your contributions will be licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

Thank you for contributing to WILDS! Your contributions help advance reproducible bioinformatics research for the entire community.