# ww-bedtools
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for using BEDTools utilities.

## Overview

This module provides reusable WDL tasks for using **BEDTools** to work with genomic intervals.

The module is designed to be a foundational component within the WILDS ecosystem, suitable for use in larger sequence analysis pipelines.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Task**: `coverage`, `intesect`, `makewindows`, `validate_outputs`
- **Workflow**: `bedtools_example` (demonstration workflow that executes all tasks)
- **Container**: `getwilds/bedtools:2.31.1`

## Tasks

### `coverage`

Calculates average read coverage over given BED intervals for a sample.

**Inputs:**

- `bed_file` (File): BED file containing genomic intervals
- `aligned_bam` (File): Input aligned and indexed BAM file
- `sample_name` (String): Name of the sample provided for output files
- `memory_gb` (Int): Memory allocation in GB (default: 16)
- `cpu_cores` (Int): Number of CPU cores (default: 2)

**Outputs:**

- `name` (String): Sample name
- `mean_coverage` (File): Output file with mean coverage values

### `intersect`

Uses BEDTools to find overlaps between the sample's BAM and a BED file.

**Inputs:**

- `bed_file` (File): BED file to intersect with
- `aligned_bam` (File): Input aligned and indexed BAM file
- `sample_name` (String): Name of the sample provided for output files
- `flags` (String): BEDTools intersect flags (default: `-header -wo`)
- `memory_gb` (Int): Memory allocation in GB (default: 16)
- `cpu_cores` (Int): Number of CPU cores (default: 2)

**Outputs:**

- `name` (String): Sample name
- `intersect_output` (File): `bedtools intersect` result file

### `makewindows`

Creates genomic windows by chromosome and counts reads overlapping each window.

**Inputs:**

- `bed_file` (File): BED file of input intervals
- `aligned_bam` (File): Input aligned BAM file
- `bam_index` (File): Index of aligned BAM file
- `reference_fasta` (File): Reference genome FASTA
- `reference_index` (File): Reference genome index
- `list_chr` (Array[String]): List of chromosomes to analyze
- `sample_name` (String): "Name of the sample provided for output files"
- `tmp_dir` (String): Temporary directory path
- `memory_gb` (Int): Memory allocation in GB (default: 24)
- `cpu_cores` (Int): Number of CPU cores (default: 10)

**Outputs:**

- `name` (String): Sample name
- `counts_bed` (File): `.tar.gz` archive of per-chromosome BED files with read counts


### `validate_outputs`

Checks existence and content of all output files by sample.

**Inputs:**

- `intersect_files` (Array[File]): Intersect output files
- `coverage_files` (Array[File]): Mean coverage files
- `window_count_files` (Array[File]): BED file tarballs with read counts
- `sample_names` (Array[String]): Sample names corresponding to each file

**Outputs:**

- `report` (File): Human-readable validation report

## Workflow: `bedtools_example`

The included example workflow processes BAM files using all BEDTools tasks and validates the results.

**Inputs:**

- `bed_file` (File): BED file of intervals
- `samples` (Array[SampleInfo]): Structs containing BAM, index, and sample name
- `reference_fasta` (File): Genome reference FASTA
- `reference_index` (File): FAI index for reference
- `intersect_flags` (String): Flags for `bedtools intersect` (default: `-header -wo`)
- `chromosomes` (Array[String]): Chromosomes for `makewindows`
- `tmp_dir` (String): Temporary path for intermediate files
- `cpus` (Int): Threads for each task (default: 8)
- `memory_gb` (Int): Memory for each task (default: 32)

**Outputs:**

- `intersect_results` (Array[File]): Output from `intersect`
- `coverage_results` (Array[File]): Output from `coverage`
- `window_count_results` (Array[File]): Output tarballs from `makewindows`
- `bedtools_validation_report` (File): Output report from `validate_outputs`

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/main/modules/ww-bedtools/ww-bedtools.wdl" as bedtools

workflow my_pipeline {
  input {
    Array[bedtools.SampleInfo] samples
    File bed_file
    File reference_fasta
    File reference_index
  }

  call bedtools.bedtools_example {
    input:
      bed_file = bed_file,
      samples = samples,
      reference_fasta = reference_fasta,
      reference_index = reference_index
  }

  output {
    File validation_report = bedtools_example.bedtools_validation_report
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-sra**: Download sequencing data prior to alignment
- **ww-bwa**: Sequence alignment

## Testing the Module

The module includes a demonstration workflow (`bedtools_example`) with support for execution on multiple WDL backends:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-bedtools.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-bedtools.wdl -i inputs.json

# Using Sprocket
sprocket run ww-bedtools.wdl inputs.json
```

### Test Input Format

```json
{
  "bedtools_example.samples": [
    {
      "name": "sample1",
      "bam": "/path/to/sample1.bam",
      "bam_index": "/path/to/sample1.bam.bai"
    }
  ],
  "bedtools_example.bed_file": "/path/to/regions.bed",
  "bedtools_example.reference_fasta": "/path/to/genome.fasta",
  "bedtools_example.reference_index": "/path/to/genome.fasta.fai",
  "bedtools_example.tmp_dir": "/tmp"
}
```

## Configuration Guidelines

### Resource Allocation

- **Memory**: 16-24 GB recommended for full human genomes; can be tuned for smaller input data (such as just one chromosome).
- **CPUs**: 2-10 cores for full human genomes. `makewindows` requires 10 cores, for the others 2 cores should be sufficient.

### Advanced Parameters

- Ensure the reference genome is pre-indexed and has an associated `.fai`
- Set `cpu_cores` and `memory_gb` based on your environment and dataset size

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker or Apptainer support
- Sufficient memory for genome indexing (varies by genome size)

## Features

- **BEDTools intersect analysis**: Screen for overlaps between two sets of genomic features
- **BEDTools coverage analysis**: Read coverage results across BED intervals
- **BEDTools makewindows output**: Per-chromosome BED files of read counts
- **Validation**: Built-in output validation and reporting
- **Modular design**: Integrates with other WILDS workflows and tools
- **Scalable**: Supports batch alignment across many samples
- **Flexible**: Customizable resource settings per task
- **Robust**: Includes error handling and reproducible output filenames

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real RNA-seq data (chromosome 22 subset for efficiency)
- Comprehensive validation of all outputs

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, contact the Fred Hutch Data Science Lab (DaSL) at [wilds@fredhutch.org](mailto:wilds@fredhutch.org) or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

To contribute to this module, please review the [WILDS Contributor Guide](https://getwilds.org/guide/) and the [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md).

## License

Distributed under the MIT License. See `LICENSE` for full details.