<table>
<tr>
  <td style="width: 400px; max-width: 400px;"><img src="https://raw.githubusercontent.com/getwilds/wilds-wdl-library/main/WILDSWDLLogo.jpeg" width="400" alt="WILDS WDL logo"></td>
  <td style="vertical-align: middle;">
    <h1>WILDS WDL Library Documentation</h1>
    <p>Technical documentation for a centralized collection of bioinformatics WDL infrastructure providing reusable, well-tested components for genomics research.</p>
  </td>
</tr>
</table>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![WDL Executors](https://img.shields.io/badge/WDL-Cromwell%20%7C%20miniWDL%20%7C%20Sprocket-blue.svg)](https://github.com/getwilds/wilds-wdl-library)
[![WDL](https://img.shields.io/badge/WDL-1.0-orange.svg)](https://openwdl.org/)

---

## Welcome to the WILDS WDL Library Documentation

This site provides comprehensive technical documentation for all WDL modules and pipelines in the WILDS WDL Library. Use the **sidebar navigation** to explore available components and view detailed task descriptions, inputs/outputs, and usage examples.

### Navigating This Documentation

**Modules** - Foundational building blocks
- Tool-specific collections of reusable WDL tasks
- Each module page shows available tasks with complete parameter documentation
- Includes test workflows demonstrating basic usage

**Pipelines** - Complete analysis workflows
- Workflows combining multiple modules into functional analysis pipelines
- Range from basic educational examples (2-3 modules) to advanced production pipelines (10+ modules)
- Complexity levels documented in each pipeline's README
- Serve as templates for custom workflows or production analyses

---

## Quick Start: Using Components

### Importing Modules into Your Workflows

All components can be imported directly via GitHub URLs:

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sra/ww-sra.wdl" as sra_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-star/ww-star.wdl" as star_tasks

workflow my_analysis {
  call sra_tasks.fastqdump { input: sra_id = "SRR12345678" }
  call star_tasks.star_align_two_pass {
    input: sample_data = { "name": "sample1", "r1": fastqdump.r1_end, "r2": fastqdump.r2_end }
  }
}
```

**Benefits of GitHub imports:**
- No local cloning required - use modules directly
- Pin to specific commits or tags for reproducibility
- Easy version switching by changing the URL
- Import only the modules you need

### Running Pipelines Directly (No Clone Required)

Thanks to GitHub URL imports, you can download and run any pipeline without cloning the entire repository:

```bash
# Download a pipeline and its example inputs
# Option 1: Use curl from the command line
curl -O https://raw.githubusercontent.com/getwilds/wilds-wdl-library/main/pipelines/ww-sra-star/ww-sra-star.wdl
curl -O https://raw.githubusercontent.com/getwilds/wilds-wdl-library/main/pipelines/ww-sra-star/inputs.json
# Option 2: Download directly from GitHub by navigating to the file and clicking the download button

# Modify inputs.json as necessary for your data, then run via the command line or PROOF's point-and-click interface
sprocket run ww-sra-star.wdl inputs.json
```

### Running Components Locally

If you prefer to clone the full repository:

```bash
# Clone the repository
git clone https://github.com/getwilds/wilds-wdl-library.git
cd wilds-wdl-library

# Run a module test workflow (no inputs needed)
cd modules/ww-star
sprocket run testrun.wdl

# Run a pipeline (modify inputs.json as necessary)
cd ../../pipelines/ww-sra-star
sprocket run ww-sra-star.wdl inputs.json
```

---

## Library Architecture

The WILDS WDL Library is organized into two complementary tiers:

### Modules (`modules/`)
**Purpose**: Foundational building blocks for larger workflows
**Content**: Individual bioinformatics tools (STAR, BWA, GATK, etc.)
**Testing**: Unit tests ensure each task functions correctly
**Usage**: Import tasks into custom workflows or run demonstration workflows

### Pipelines (`pipelines/`)
**Purpose**: Functional pipelines ranging from educational examples to production-ready analyses
**Content**: Multiple modules combined into analysis workflows of varying complexity
**Complexity Levels**: Basic (2-3 modules), Intermediate (4-6 modules), Advanced (10+ modules)
**Testing**: Integration tests verify modules work together seamlessly
**Usage**: Templates for common workflows, learning examples, or production analyses

---

## Supported WDL Executors

All components are tested with multiple WDL executors to ensure broad compatibility:

- **[Cromwell](https://cromwell.readthedocs.io/)** - Production-grade workflow engine
- **[miniWDL](https://github.com/chanzuckerberg/miniwdl)** - Lightweight local execution
- **[Sprocket](https://sprocket.bio/)** - Modern WDL executor with enhanced features

---

## Container Images

All tasks use versioned, tested Docker images from the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library), ensuring reproducible execution across different computing environments.

---

## Getting Help

- **Documentation Issues**: Found something unclear or incorrect? [Report an issue](https://github.com/getwilds/wilds-wdl-library/issues)
- **General Questions**: Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- **Additional Resources**:
  - [WILDS Guide](https://getwilds.org/guide/) - Comprehensive guides and best practices
  - [Contributing Guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) - How to contribute
  - [GitHub Repository](https://github.com/getwilds/wilds-wdl-library) - Source code and development

---

**Ready to explore?** Use the sidebar to browse available modules and pipelines. Each component page provides complete technical documentation including task signatures, parameter descriptions, and usage examples.
