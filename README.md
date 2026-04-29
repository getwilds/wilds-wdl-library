<table>
<tr>
  <td style="width: 400px; max-width: 400px;"><img src="https://raw.githubusercontent.com/getwilds/wilds-wdl-library/main/WILDSWDLLogo.jpeg" width="400" alt="WILDS WDL logo"></td>
  <td style="vertical-align: middle;">
    <h1>WILDS WDL Library</h1>
    <p>A centralized collection of bioinformatics WDL infrastructure providing reusable, well-tested components that can be combined to create powerful computational pipelines for genomics research.</p>
  </td>
</tr>
</table>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Project Status: Stable – Useable, full support, open to feedback, stable API.](https://getwilds.org/badges/badges/stable.svg)](https://getwilds.org/badges/#stable)
[![WDL Executors](https://img.shields.io/badge/WDL-Cromwell%20%7C%20miniWDL%20%7C%20Sprocket-blue.svg)](https://github.com/getwilds/wilds-wdl-library)
[![GitHub Release](https://img.shields.io/github/v/release/getwilds/wilds-wdl-library)](https://github.com/getwilds/wilds-wdl-library/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19339693.svg)](https://doi.org/10.5281/zenodo.19339693)<br>
[![Module Tests](https://github.com/getwilds/wilds-wdl-library/actions/workflows/modules-testrun.yml/badge.svg)](https://github.com/getwilds/wilds-wdl-library/actions/workflows/modules-testrun.yml)
[![Pipeline Tests](https://github.com/getwilds/wilds-wdl-library/actions/workflows/pipelines-testrun.yml/badge.svg)](https://github.com/getwilds/wilds-wdl-library/actions/workflows/pipelines-testrun.yml)
[![Linting](https://github.com/getwilds/wilds-wdl-library/actions/workflows/linting.yml/badge.svg)](https://github.com/getwilds/wilds-wdl-library/actions/workflows/linting.yml)

## Overview

The WILDS WDL Library consolidates bioinformatics workflows into a single, well-organized repository that serves as both a collection of production-ready tools and a demonstration of WDL best practices. Rather than maintaining separate repositories for each workflow, this library promotes modularity and reusability through a two-tier architecture.

## Library Architecture

The library is organized into two complementary levels:

### **Modules** (`modules/`)
Tool-specific collections of reusable WDL tasks with comprehensive testing.
- **Purpose**: Foundational building blocks for larger workflows
- **Content**: Individual bioinformatics tools (STAR, BWA, GATK, etc.)
- **Testing**: Unit tests ensure each task functions correctly over time
- **Usage**: Import tasks into custom workflows or run demonstration workflows

### **Pipelines** (`pipelines/`)
Complete analysis workflows combining multiple modules.
- **Purpose**: Functional pipelines ranging from educational examples to production-ready analyses
- **Content**: Multiple modules combined into analysis workflows of varying complexity
- **Complexity Levels**: Basic (2-3 modules), Intermediate (4-6 modules), Advanced (10+ modules)
- **Testing**: Integration tests verify modules work together seamlessly
- **Usage**: Templates for common workflows, learning examples, or production analyses

## Quick Start

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

### Using the Full Repository

If you want to explore multiple components or contribute:

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

### Importing into Your Workflows

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

WILDS pipelines use GitHub URLs for imports, providing several advantages:

- **No local cloning required**: Use modules directly without downloading the repository
- **Version control**: Pin to specific commits or tags for reproducibility
- **Easy updates**: Switch between versions by changing the URL
- **Modular usage**: Import only the modules you need

## Supported Executors

All components are tested with multiple WDL executors:
- **[Sprocket](https://sprocket.bio/)**: Modern WDL executor with enhanced features
- **[Cromwell](https://cromwell.readthedocs.io/)**: Production-grade workflow engine
- **[miniWDL](https://github.com/chanzuckerberg/miniwdl)**: Lightweight local execution

## For Fred Hutch Users

Fred Hutch researchers can use [PROOF](https://sciwiki.fredhutch.org/datademos/proof-how-to/) to submit workflows directly to the on-premise HPC cluster. This provides a user-friendly interface for researchers unfamiliar with command-line tools while leveraging the power of the institutional computing resources.

**Cromwell Configuration**: PROOF users can customize workflow execution using Cromwell options. See [cromwell-options.json](cromwell-options.json) for example configurations including call caching, output directories, and more. For detailed information, refer to the [Cromwell workflow options documentation](https://cromwell.readthedocs.io/en/stable/wf_options/Overview/).

**Platform-Specific Configurations**: Some pipelines include optional platform-specific configurations (e.g., `.cirro/` directories) for execution on cloud platforms like [Cirro](https://cirro.bio/). These configurations are self-contained within each pipeline directory.

## Quality Assurance

### Automated Testing
- **Continuous Integration**: All components tested on every pull request
- **Multi-Executor Validation**: Ensures compatibility across different WDL engines
- **Real Data Testing**: Uses authentic bioinformatics datasets for validation
- **Scheduled Monitoring**: Weekly checks detect infrastructure changes
- **HPC Monthly Test Runs**: Resource-intensive modules excluded from GitHub CI are validated monthly on the Fred Hutch HPC, with results posted to a [tracking issue](https://github.com/getwilds/wilds-wdl-library/issues) for visibility

### Standards and Best Practices
- **Standardized Structure**: Consistent organization across all components
- **Container Management**: Versioned, tested Docker images from the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)
- **Documentation Standards**: Comprehensive README files and inline documentation
- **Version Control**: Semantic versioning and careful dependency management

## AI Disclosure

Large language models (primarily [Claude](https://www.anthropic.com/claude)) have been used to assist with prototyping new modules and pipelines in this repository, including initial drafts of task scaffolding, command structures, and test workflows. All AI-generated code is reviewed, validated, and modified as necessary by human researchers through linting tools, multi-executor test runs (Cromwell, miniWDL, Sprocket), and peer code review. No architectural or scientific decisions are made by AI tools.

## Contributing

We welcome contributions at all levels:

### Adding New Modules
1. Focus on high-utility bioinformatics tools
2. Follow the standard module structure
3. Include comprehensive tests and validation
4. Provide detailed documentation

### Creating Pipelines
1. Combine existing modules (prefer existing modules over new tasks)
2. Demonstrate common analysis patterns
3. Include realistic test datasets
4. Document complexity level and integration approaches

### Improving Documentation
- Enhance existing README files
- Add usage examples and tutorials
- Improve inline code documentation
- Contribute to the [WILDS documentation site](https://getwilds.org/)

See our [Contributing Guidelines](.github/CONTRIBUTING.md) for detailed information.

## Development Roadmap

### Current Focus
- Growing the library through the **WILDS WDL Development Program**, a collaborative initiative where Fred Hutch researchers partner with the OCDO Research Informatics team to convert their computational analyses into production-quality WDL workflows
- Expanding module and pipeline coverage based on researcher needs identified through the development program
- Enhancing cloud platform integration (e.g., [Cirro](https://cirro.bio/)) for broader execution options

### Future Plans
- Continued researcher-driven development of modules and pipelines through the WILDS WDL Development Program
- Additional advanced pipelines for publication-ready analyses
- Enhanced documentation, tutorials, and community outreach

## Support

- **Issues and Bug Reports**: [GitHub Issues](https://github.com/getwilds/wilds-wdl-library/issues)
- **General Questions**: Contact the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org
- **Documentation**: [Contributing Guidelines](.github/CONTRIBUTING.md)
- **Fred Hutch Users**: [Scientific Computing Wiki](https://sciwiki.fredhutch.org/)

## Citation

If you use the WILDS WDL Library in your research, please cite it. The DOI below always resolves to the latest release; see the [Zenodo record](https://doi.org/10.5281/zenodo.19339693) for version-specific DOIs.

> Firman, T., Bishop, E., et al. (2026). WILDS WDL Library [Computer software]. GitHub. https://doi.org/10.5281/zenodo.19339693

BibTeX:

```bibtex
@software{wilds_wdl_library,
  author    = {Firman, Taylor and Bishop, Emma and others},
  title     = {WILDS WDL Library},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.19339693},
  url       = {https://doi.org/10.5281/zenodo.19339693}
}
```

A machine-readable citation is also available in [CITATION.cff](CITATION.cff), which GitHub uses to power the "Cite this repository" button in the sidebar.

## Dockstore

All modules and pipelines in this library are published on [Dockstore](https://dockstore.org/), enabling discovery, import, and citation of individual workflows. Browse the library on Dockstore to find modules and pipelines with detailed descriptions, author information, and DOIs.

## Related Resources

- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[Fred Hutch SciWiki](https://sciwiki.fredhutch.org/)**: Institutional computing resources and tutorials

## License

Distributed under the MIT License. See `LICENSE` for details.
