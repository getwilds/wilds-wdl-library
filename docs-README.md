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

This site provides comprehensive technical documentation for all WDL modules, vignettes, and workflows in the WILDS WDL Library. Use the **sidebar navigation** to explore available components and view detailed task descriptions, inputs/outputs, and usage examples.

### Navigating This Documentation

**Modules** - Foundational building blocks
- Tool-specific collections of reusable WDL tasks
- Each module page shows available tasks with complete parameter documentation
- Includes test workflows demonstrating basic usage

**Vignettes** - Integration examples
- Compact workflows combining 2-3 modules into common analysis patterns
- Demonstrates how modules work together
- Serves as templates for custom workflows

**Workflows** - Production pipelines
- Complete, publication-ready analysis pipelines
- Complex workflows combining multiple modules
- Ready for production use with minimal customization

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

### Running Components Locally

```bash
# Clone the repository
git clone https://github.com/getwilds/wilds-wdl-library.git
cd wilds-wdl-library

# Run a module test workflow (no inputs needed)
cd modules/ww-star
sprocket run testrun.wdl

# Run a vignette (update inputs.json as needed)
cd ../../vignettes/ww-sra-star
sprocket run ww-sra-star.wdl inputs.json

# Run a full workflow (update inputs.json as needed)
cd ../../workflows/ww-leukemia
sprocket run ww-leukemia.wdl inputs.json
```

---

## Library Architecture

The WILDS WDL Library is organized into three complementary tiers:

### Modules (`modules/`)
**Purpose**: Foundational building blocks for larger workflows
**Content**: Individual bioinformatics tools (STAR, BWA, GATK, etc.)
**Testing**: Unit tests ensure each task functions correctly
**Usage**: Import tasks into custom workflows or run demonstration workflows

### Vignettes (`vignettes/`)
**Purpose**: Educational examples of module integration
**Content**: 2-3 modules combined into standard analysis patterns
**Testing**: Integration tests verify modules work together
**Usage**: Templates for common workflows or learning examples

### Workflows (`workflows/`)
**Purpose**: End-to-end analyses suitable for research publications
**Content**: Complex workflows combining multiple modules and custom logic
**Testing**: Comprehensive validation with realistic datasets
**Usage**: Production analyses requiring minimal customization

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

**Ready to explore?** Use the sidebar to browse available modules, vignettes, and workflows. Each component page provides complete technical documentation including task signatures, parameter descriptions, and usage examples.
