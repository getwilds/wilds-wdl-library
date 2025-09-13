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
- [Pull Request Process](#pull-request-process)
- [Code of Conduct](#code-of-conduct)

## Getting Started

Before contributing code changes, please:

1. **[Fork the repository](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo)** to your GitHub account

2. **Set up your development environment** with the required tools:
   - For local testing:
       - [miniWDL](https://miniwdl.readthedocs.io/en/latest/getting_started.html#install-miniwdl)
       - [sprocket](https://sprocket.bio/installation.html)
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

1. **`ww-toolname.wdl`** - Main WDL workflow file named for the tool it uses
2. **`README.md`** - Comprehensive documentation


**Your WDL file must include:**

- **Version declaration**: Use WDL version 1.0
- **Example workflow**: A `toolname_example` workflow that uses all tasks for testing
- **Task definitions**: Individual tasks with proper resource requirements
- **Proper imports**: Use GitHub URLs for dependencies on existing WILDS WDL modules (e.g. for downloading test data)
- **Metadata documentation**: Describe properties of the workflow and tasks (e.g. inputs, outputs)

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

## Workflow Development Guidelines

**Workflows represent the highest tier and must meet production standards:**

- **Comprehensive validation**: Extensive testing with realistic datasets
- **Scalability**: Support for large datasets and batch processing
- **Documentation**: Thorough documentation to aid in publication

## Testing Requirements

### Local Tests

Make sure you have these installed:

- [miniWDL](https://miniwdl.readthedocs.io/en/latest/getting_started.html#install-miniwdl)
- [sprocket](https://sprocket.bio/installation.html)
- [Docker Desktop](https://www.docker.com/get-started/) for container execution

Test your WDL locally before submitting a PR:

```bash
cd modules/ww-toolname

# Linting with miniwdl
miniwdl check ww-toolname.wdl

# Linting with sprocket (ignoring things we don't care about)
sprocket lint \
  -e TodoComment \
  -e ContainerUri \
  -e TrailingComma \
  -e CommentWhitespace \
  -e UnusedInput \
  ww-toolname.wdl

# Test running
miniwdl run ww-toolname.wdl
sprocket run ww-toolname.wdl
```

For automated local testing be sure to have [uv](https://docs.astral.sh/uv/getting-started/installation/) installed also and use our Makefile.

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