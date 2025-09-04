# Contributing to WILDS WDL Library

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

1. **Fork the repository** to your GitHub account
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/your-username/wilds-wdl-library.git
   cd wilds-wdl-library
   ```
3. **Set up development environment** with required tools:
   - [miniWDL](https://github.com/chanzuckerberg/miniwdl) for local testing
   - [Docker](https://www.docker.com/) for container execution
   - A WDL-compatible workflow engine ([Cromwell](https://cromwell.readthedocs.io/), [miniWDL](https://github.com/chanzuckerberg/miniwdl), or [Sprocket](https://sprocket.bio/))

4. **Create a feature branch** for your contribution:
   ```bash
   git checkout -b new-branch-name
   ```

   Commit and push your changes to this new branch.

5. **Create a pull request (PR)** to merge your contribution branch into `main`.
    *   The title of your PR should briefly describe the change.
    *   If your contribution resolves an issue, the body of your PR should contain `Fixes #issue-number`.

## Repository Structure

The WILDS WDL Library follows a three-tier architecture:

```
wilds-wdl-library/
├── modules/           # Individual tool modules (Tier 1)
│   └── ww-toolname/
│       ├── ww-toolname.wdl
│       ├── options.json
│       ├── inputs.json
│       └── README.md
├── vignettes/         # Educational workflow examples (Tier 2)
│   └── ww-example-name/
│       ├── ww-example-name.wdl
│       ├── inputs.json
│       └── README.md
├── workflows/         # Production-ready pipelines (Tier 3)
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
- Add usage tutorials and best practices
- Contribute to the [WILDS documentation site](https://getwilds.org/)

### 3. Module Contributions (Tier 1)

- Focus on one high-utility bioinformatics tool
- Follow standardized module structure
- Include comprehensive testing and validation

### 4. Vignette Contributions (Tier 2)

- Combine 2-3 existing modules
- Demonstrate common analysis patterns
- Provide educational value for the community

### 5. Workflow Contributions (Tier 3)

- Create publication-ready analysis pipelines
- Integrate multiple modules with custom logic
- Include comprehensive validation with realistic datasets

## Module Development Guidelines

### Required Files

Each module must contain exactly these four files:

1. **`ww-toolname.wdl`** - Main WDL workflow file named for the tool it uses
2. **`options.json`** - Runtime configuration options
3. **`inputs.json`** - Example input parameters
4. **`README.md`** - Comprehensive documentation

### WDL File Requirements

Your WDL file must include:

- **Version declaration**: Use WDL version 1.0 or higher
- **Task definitions**: Individual tool tasks with proper resource requirements
- **Test workflow**: A `test_data` workflow that exercises all tasks
- **Validation task**: A `validate_outputs` task that produces a text file report
- **Proper imports**: Use GitHub URLs for any external dependencies on existing WILDS WDL modules

Example structure:
```wdl
version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/main/modules/ww-testdata/ww-testdata.wdl" as testdata

workflow test_data {
  call testdata.get_test_data
  call your_main_task { input: data = get_test_data.output }
  call validate_outputs { input: results = your_main_task.output }
}

task your_main_task {
  # Your tool implementation
}

task validate_outputs {
  input {
    # Input validation parameters
  }
  command <<<
    # Generate validation report
    echo "Validation report for module outputs" > validation_report.txt
    # Add specific validation logic
  >>>
  output {
    File validation_report = "validation_report.txt"
  }
}
```

### Container Requirements

- Use containers from the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library) when available
- If creating new containers, follow WILDS container standards and consider contributing to the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library).
- Specify exact container versions (avoid `latest` tags)
- Document container dependencies in the README

### Input/Output Specifications

- Use descriptive parameter names
- Provide comprehensive documentation for all inputs and outputs
- Include optional parameters with sensible defaults
- Support both single samples and batch processing where applicable

## Vignette Development Guidelines

### Purpose and Scope

Vignettes should:
- Combine 2-3 existing modules from the library
- Demonstrate realistic analysis workflows
- Serve as educational templates
- Use publicly available test data

### Required Components

- **Integration focus**: Show how modules work together seamlessly
- **Educational value**: Include explanatory comments and documentation
- **Realistic examples**: Use authentic bioinformatics datasets
- **Clear documentation**: Explain the biological/analytical context

### No New Tasks Rule

Vignettes should **only** combine existing modules - do not create new task definitions. If you need new functionality, contribute it as a module first.

## Workflow Development Guidelines

### Production-Ready Standards

Workflows represent the highest tier and must meet production standards:

- **Comprehensive validation**: Extensive testing with realistic datasets
- **Scalability**: Support for large datasets and batch processing
- **Documentation**: Thorough documentation

### Publication Readiness

- Include detailed documentation suitable that would be useful for a methods section
- Provide and resource requirements and performance benchmarks if possible
- Document software versions and parameter choices

## Testing Requirements

### Automated Testing

All contributions must pass our automated testing pipeline which executes on a PR via GitHub Actions:

- **Multi-executor validation**: Tested with Cromwell, miniWDL, and Sprocket
- **Container verification**: All Docker images must be accessible and functional
- **Syntax validation**: WDL syntax and structure validation
- **Integration testing**: Cross-module compatibility testing

### Local Testing

Before submitting, test your contribution locally:

```bash
# Test with miniWDL
cd modules/your-module
miniwdl run ww-your-module.wdl -i inputs.json
```

### Test Data Requirements

- Use the `ww-testdata` module for standardized test datasets
- If you need additional test datasets, modify the `ww-testdata` module also
- Include small, representative test files in your examples

## Documentation Standards

### README Requirements

Each contribution must include a comprehensive README.md:

#### Module README Template

````markdown
# ww-toolname

Brief description of the tool and its purpose.

## Overview

Detailed description of what the module does, its biological context, and use cases.

## Inputs

| Parameter | Type | Description | Required |
|-----------|------|-------------|----------|
| input_file | File | Description of input | Yes |
| parameter | String | Description of parameter | No (default: value) |

## Outputs

| Output | Type | Description |
|--------|------|-------------|
| output_file | File | Description of output |

## Usage

### Basic Usage
```bash
miniwdl run ww-toolname.wdl -i inputs.json
```

### In Other Workflows
```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/main/modules/ww-toolname/ww-toolname.wdl" as toolname
```

## Container

- **Image**: `getwilds/toolname:version`
- **Source**: [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)

## Testing

Describe the test workflow and validation approach.

## Citation

Include appropriate citations for the underlying tool.
````

### Code Documentation

- Use clear, descriptive variable names
- Include inline comments for complex logic
- Document runtime requirements and resource usage
- Explain parameter choices and defaults

## Pull Request Process

### Before Submitting Code Changes

1. **Test thoroughly**: Ensure all tests pass locally
2. **Update documentation**: Including the repository README as appropriate (e.g. when adding a new module)
3. **Follow naming conventions**: Use the `ww-` prefix for all WDL files
4. **Check file structure**: Ensure all required files are present (e.g. `inputs.json` when submitting a new module)

### Submission Process

1. **Create descriptive PR title**:
   - Examples: `Add BWA alignment module`, `RNA-seq analysis example vingette`

2. **Fill out PR template**: Provide detailed information about your contribution

3. **Link related issues**: Reference any GitHub issues your PR addresses

4. **Request reviews**: Tag Emma Bishop (@emjbishop) or Taylor Firman (@tfirman)

### Review Criteria

Your PR will be evaluated on:

- **Functionality**: Does it work as intended?
- **Testing**: Are tests comprehensive and passing?
- **Documentation**: Is documentation clear and complete?
- **Standards compliance**: Does it follow WILDS conventions?
- **Code quality**: Is the WDL code well-structured and readable?

### CI/CD Validation

All PRs must pass automated checks:

- WDL syntax validation
- Container accessibility testing
- Multi-executor compatibility testing
- Documentation completeness checks
- Security scanning

## Development Best Practices

### WDL Best Practices

- **Use appropriate WDL version**: Prefer version 1.0 or newer
- **Resource specification**: Always specify CPU, memory, and disk requirements
- **Error handling**: Include meaningful error messages and exit codes
- **Modularity**: Create reusable tasks that can be imported by other workflows
- **Parameter validation**: Validate inputs when possible

### Container Best Practices

- **Version pinning**: Always specify exact container versions
- **Minimal images**: Use lightweight base images when possible
- **Security**: Scan containers for vulnerabilities
- **Documentation**: Document all installed tools and versions

### Performance Considerations

- **Resource estimation**: Provide accurate resource requirements
- **Scalability**: Design for both small test runs and large production datasets
- **Efficiency**: Optimize for computational and storage efficiency
- **Monitoring**: Include logging for debugging and monitoring

## Getting Help

### Communication Channels

- **GitHub Issues**: For bug reports and feature requests
- **Email**: Contact the Fred Hutch Data Science Lab at [wilds@fredhutch.org](mailto:wilds@fredhutch.org)
- **Documentation**: [WILDS Guide](https://getwilds.org/guide/)
- **Fred Hutch Users**: [Scientific Computing Wiki](https://sciwiki.fredhutch.org/)

### Mentorship

New contributors are welcome! If you're new to WDL or bioinformatics workflows:

- Review our [WDL 101 course materials](https://github.com/getwilds/wdl-101)
- Check out existing modules for examples
- Don't hesitate to ask questions in issues or via email. If you have a `uw.edu` or `fredhutch.org` email you can also ask questions in our `fh-data` [slack workspace]([https://hutchdatascience.org/joinslack/)
- Consider starting with documentation contributions

## Recognition

Contributors will be acknowledged in:

- Repository contributor lists
- Module-specific acknowledgments
- WILDS project publications where appropriate
- Community showcases and presentations

## Code of Conduct

By participating in this project, you agree to abide by our code of conduct:

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