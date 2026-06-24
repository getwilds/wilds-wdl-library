# ww-mosdepth Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for calculating sequencing coverage depth using mosdepth.

## Overview

This module wraps `mosdepth`, a fast tool for calculating sequencing coverage depth from BAM or CRAM alignments at per-base, per-region, or per-window resolution. It provides a standardized way to generate coverage statistics and depth files.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-mosdepth.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `calculate_depth`

Calculates sequencing coverage depth from BAM or CRAM alignments.

**Inputs:**
- `sample_name` (String): Name identifier for the sample
- `input_bam` (File): Input BAM or CRAM file
- `input_bam_index` (File): Index file for the input BAM/CRAM
- `ref_fasta` (File, optional): Reference genome FASTA file (required for CRAM inputs)
- `regions_bed` (File, optional): BED file specifying regions of interest
- `window_size` (Int, default=100): Window size for windowed depth output (in bp)
- `cpu_cores` (Int, default=2): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=4): Memory allocated for the task in GB

**Outputs:**
- `depth_per_base` (File): Per-base coverage depth in BED format
- `depth_summary` (File): Summary of coverage depth statistics
- `region_depth` (File): Coverage depth for specified regions or windows

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-mosdepth/ww-mosdepth.wdl" as mosdepth_tasks

workflow my_coverage_analysis {
  input {
    File input_bam
    File input_bai
    File ref_fasta
  }
  
  call mosdepth_tasks.calculate_depth {
    input: {
      sample_name: "sample_1",
      input_bam: input_bam,
      input_bam_index: input_bai,
      ref_fasta: ref_fasta
    }
  }
}
```

### Integration Examples

This module pairs well with other WILDS modules:
- **ww-testdata**: Automatic provisioning of test BAM and reference data
- **ww-bwa** / **ww-bowtie2** / **ww-star**: For producing the aligned BAM/CRAM inputs
- **ww-samtools**: For preprocessing (sorting, indexing, filtering) prior to depth calculation

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl
```

## Docker Container

This module uses the `getwilds/mosdepth:0.3.14` container image, which includes:
- The `mosdepth` binary (v0.3.14)
- All necessary system dependencies for coverage depth calculation

## Citation

If you use mosdepth in your research, please cite the original authors:

> Pedersen BS, Quinlan AR. Mosdepth: quick coverage depth calculation for genomes and exomes. Bioinformatics. 2018 Mar 1;34(5):867-868.
> DOI: 10.1093/bioinformatics/btx699

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 2 cores
- **Memory**: 4 GB

### Resource Scaling
- `cpu_cores`: Increase for larger BAM files or more chromosomes
- `memory_gb`: Increase based on the size of the reference genome and the input BAM

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Sufficient computational resources (mosdepth is generally lightweight, but scales with BAM size)

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data (chromosome 1 subset for efficiency)
- Validation of expected output files

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to mosdepth usage or configuration, please refer to the [mosdepth repository](https://github.com/brentp/mosdepth).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
