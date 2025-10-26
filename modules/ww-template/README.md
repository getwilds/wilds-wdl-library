# ww-template Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL template module demonstrating best practices for creating new modules. This functional template performs simple "hello world" operations purely for testing purposes and serves as a copy-paste starting point for new tool integrations.

## Overview

This module provides a complete template for creating new WILDS WDL modules. It includes all the standard components and patterns used across the WILDS ecosystem, with simple echo-based functionality that demonstrates key concepts like testing and integration.

**For Contributors**: Use this module as a starting point for wrapping new bioinformatics tools. Simply copy the structure and replace the echo commands with your tool's specific functionality.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-template.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `process_sample`

Simple template processing task that creates a hello world output file.

**Inputs:**
- `sample_name` (String): Name identifier for the sample
- `input_file` (File): Input file (any file type works for this template)
- `cpu_cores` (Int, default=1): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=4): Memory allocated for the task in GB

**Outputs:**
- `output_file` (File): Simple text file with hello world message and sample information


## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-template/ww-template.wdl" as template_tasks

struct TemplateSample {
    String name
    File input_file
}

workflow my_analysis_pipeline {
  input {
    Array[TemplateSample] samples
  }
  
  scatter (sample in samples) {
    call template_tasks.process_sample {
      input:
        sample_name = sample.name,
        input_file = sample.input_file
    }
  }
  
  output {
    Array[File] output_files = process_sample.output_file
  }
}
```

### Advanced Usage Examples

**Custom resource allocation:**
```wdl
call template_tasks.process_sample {
  input:
    sample_name = "large_sample",
    input_file = large_input_file,
    cpu_cores = 4,
    memory_gb = 16
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-testdata**: Automatic provisioning of test data for demonstrations
- **Other WILDS modules**: Can be used as a preprocessing step or combined with analysis modules

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint template_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads test FASTQ data using `ww-testdata`
2. Processes the test data with the simple hello world functionality
3. Demonstrates the module's tasks in a realistic workflow context

## Docker Container

This module uses the `getwilds/bwa:0.7.17` container image, which includes:
- BWA aligner (not used in template, just demonstrates container usage)
- Basic Unix text processing tools
- All necessary system dependencies for demonstration purposes

## Citation

If your module wraps a published tool or uses published reference data, include citation information here to ensure proper attribution.

**Example format:**

> Tool Name: [Paper Title]
> Authors et al. (Year)
> Journal Name, Volume(Issue), pages
> DOI: [doi link]

**Template Citation:**
This template itself doesn't require citation, but the specific tools you wrap should be properly cited according to their authors' preferences.

## Template Information

- **Purpose**: Simple template for creating new WILDS WDL modules
- **Functionality**: Basic "hello world" output generation
- **Demo Data**: Uses test FASTQ files from ww-testdata module

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 1 core
- **Memory**: 4 GB
- **Runtime**: Less than 1 minute per sample for demo data

### Resource Scaling
This template uses minimal resources by design:
- `cpu_cores`: Can be increased for CPU-intensive tools
- `memory_gb`: Can be increased for memory-intensive tools
- Resources should be adjusted based on your specific tool's requirements


## Creating Your Own Module

To use this template for a new tool:

1. **Copy the template structure**:
   ```bash
   cp -r modules/ww-template modules/ww-yourtool
   ```

2. **Update the main WDL file (`ww-yourtool.wdl`)**:
   - Rename `ww-template.wdl` to `ww-yourtool.wdl`
   - Update task definitions with your tool's specific functionality
   - Replace the `echo` commands with your tool's actual commands
   - Update Docker image to your tool's container
   - Add any additional output files your tool generates

3. **Update the test workflow (`testrun.wdl`)**:
   - Update imports to reference your new module (`ww-yourtool.wdl`)
   - Rename the workflow from `template_example` to `yourtool_example`
   - Modify struct definitions for your tool's specific inputs if needed
   - Update the test workflow to call all your module's tasks
   - Consider adding a validation task to verify output correctness

4. **Update documentation**:
   - Customize README.md with your tool's information
   - Update meta descriptions and parameter documentation in the WDL files
   - Add tool-specific usage examples
   - Add citation information if applicable

5. **Test thoroughly**:
   - Run the test workflow (`testrun.wdl`) to ensure functionality
   - Test with real data for your use case
   - Verify both linting and execution tests pass

## Contributing

To improve this template or report issues:
1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Make your changes following WILDS conventions
3. Test thoroughly with the demonstration workflow
4. Submit a pull request with detailed documentation

## Support and Feedback

For questions about this template or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines

## Related Resources

- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
