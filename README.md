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
[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![WDL Executors](https://img.shields.io/badge/WDL-Cromwell%20%7C%20miniWDL%20%7C%20Sprocket-blue.svg)](https://github.com/getwilds/wilds-wdl-library)
[![WDL](https://img.shields.io/badge/WDL-1.0-orange.svg)](https://openwdl.org/)<br>
[![Module Tests](https://github.com/getwilds/wilds-wdl-library/actions/workflows/modules-testrun.yml/badge.svg)](https://github.com/getwilds/wilds-wdl-library/actions/workflows/modules-testrun.yml)
[![Vignette Tests](https://github.com/getwilds/wilds-wdl-library/actions/workflows/vignettes-testrun.yml/badge.svg)](https://github.com/getwilds/wilds-wdl-library/actions/workflows/vignettes-testrun.yml)
[![Linting](https://github.com/getwilds/wilds-wdl-library/actions/workflows/linting.yml/badge.svg)](https://github.com/getwilds/wilds-wdl-library/actions/workflows/linting.yml)

## Overview

The WILDS WDL Library consolidates bioinformatics workflows into a single, well-organized repository that serves as both a collection of production-ready tools and a demonstration of WDL best practices. Rather than maintaining separate repositories for each workflow, this library promotes modularity and reusability through a three-tier architecture.

## Library Architecture

The library is organized into three complementary levels:

### **Modules** (`modules/`)
Tool-specific collections of reusable WDL tasks with comprehensive testing.
- **Purpose**: Foundational building blocks for larger workflows
- **Content**: Individual bioinformatics tools (STAR, BWA, GATK, etc.)
- **Testing**: Unit tests ensure each task functions correctly over time
- **Usage**: Import tasks into custom workflows or run demonstration workflows

### **Vignettes** (`vignettes/`)
Compact workflows demonstrating common bioinformatics patterns.
- **Purpose**: Educational examples of module integration
- **Content**: 2-3 modules combined into standard analysis patterns
- **Testing**: Integration tests verify modules work together seamlessly
- **Usage**: Templates for common workflows or learning examples

### **Workflows** (`workflows/`)
Complete, publication-ready analysis pipelines.
- **Purpose**: End-to-end analyses suitable for research publications
- **Content**: Complex workflows combining multiple modules and custom logic
- **Testing**: Comprehensive validation with realistic datasets
- **Usage**: Production analyses requiring minimal customization

## Quick Start

### Using Existing Components

```bash
# Clone the repository
git clone https://github.com/getwilds/wilds-wdl-library.git
cd wilds-wdl-library

# Run a module test workflow (no inputs needed)
cd modules/ww-star
miniwdl run ww-star.wdl

# Run a complete vignette (update inputs json as needed)
cd ../../vignettes/ww-sra-star
miniwdl run ww-sra-star.wdl -i inputs.json

# Run a full workflow (update inputs json as needed)
cd ../../workflows/ww-leukemia
miniwdl run ww-leukemia.wdl -i inputs.json
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

WILDS vignettes and workflows use GitHub URLs for imports, providing several advantages:

- **No local cloning required**: Use modules directly without downloading the repository
- **Version control**: Pin to specific commits or tags for reproducibility
- **Easy updates**: Switch between versions by changing the URL
- **Modular usage**: Import only the modules you need

## Supported Executors

All components are tested with multiple WDL executors:
- **[Cromwell](https://cromwell.readthedocs.io/)**: Production-grade workflow engine
- **[miniWDL](https://github.com/chanzuckerberg/miniwdl)**: Lightweight local execution
- **[Sprocket](https://sprocket.bio/)**: Modern WDL executor with enhanced features

## For Fred Hutch Users

Fred Hutch researchers can use [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit workflows directly to the on-premise HPC cluster. This provides a user-friendly interface for researchers unfamiliar with command-line tools while leveraging the power of the institutional computing resources.

## Quality Assurance

### Automated Testing
- **Continuous Integration**: All components tested on every pull request
- **Multi-Executor Validation**: Ensures compatibility across different WDL engines
- **Real Data Testing**: Uses authentic bioinformatics datasets for validation
- **Scheduled Monitoring**: Weekly checks detect infrastructure changes

### Standards and Best Practices
- **Standardized Structure**: Consistent organization across all components
- **Container Management**: Versioned, tested Docker images from the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)
- **Documentation Standards**: Comprehensive README files and inline documentation
- **Version Control**: Semantic versioning and careful dependency management

## Contributing

We welcome contributions at all levels:

### Adding New Modules
1. Focus on high-utility bioinformatics tools
2. Follow the standard module structure
3. Include comprehensive tests and validation
4. Provide detailed documentation

### Creating Vignettes
1. Combine existing modules (no new tasks)
2. Demonstrate common analysis patterns
3. Include realistic test datasets
4. Document integration approaches

### Improving Documentation
- Enhance existing README files
- Add usage examples and tutorials
- Improve inline code documentation
- Contribute to the [WILDS documentation site](https://getwilds.org/)

See our [Contributing Guidelines](.github/CONTRIBUTING.md) and the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed information.

## Development Roadmap

### Current Focus
- Expanding the module collection with high-priority tools (GATK variant calling, additional alignment tools)
- Converting existing WILDS workflows to the modular architecture
- Enhancing testing infrastructure and validation

### Future Plans
- Complete workflow tier with publication-ready pipelines
- Enhanced integration with Fred Hutch computing infrastructure
- Community-contributed modules and workflows
- Advanced documentation and tutorial content

## Support

- **Issues and Bug Reports**: [GitHub Issues](https://github.com/getwilds/wilds-wdl-library/issues)
- **General Questions**: Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- **Documentation**: [WILDS Guide](https://getwilds.org/guide/)
- **Fred Hutch Users**: [Scientific Computing Wiki](https://sciwiki.fredhutch.org/)

## Related Resources

- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[Fred Hutch SciWiki](https://sciwiki.fredhutch.org/)**: Institutional computing resources and tutorials

## License

Distributed under the MIT License. See `LICENSE` for details.
