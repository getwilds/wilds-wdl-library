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

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl
```

## Docker Container

This module uses the `getwilds/mosdepth:0.3.14` container image, which provides the `mosdepth` binary and necessary runtime dependencies.

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
