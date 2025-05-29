# ww-gatk
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for sequence alignment using the Burrows-Wheler Alignment tool (BWA).

## Overview

This module provides reusable WDL tasks for aligning reads to a reference using BWA-MEM.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `bwa_mem`
- **Workflow**: `bwa_example` (demonstration workflow that executes all tasks)
- **Container**: `getwilds/bwa:0.7.17`

## Tasks

### `bwa_mem`
Calls germline variants.

**Inputs:**
- `sample_data` (SampleInfo): Sample information struct (name, R1 and R2 FASTQs)
- `reference_fasta` (File): Reference genome FASTA file
- `cpu_cores` (Int): Total number of CPU cores allocated for the task
- `memory_gb` (Int): Memory allocated for the task in GB

**Outputs:**
- `sorted_bam` (File): Sorted BWA-MEM alignment output BAM file

## Features

The module is designed to be a foundational component within the WILDS ecosystem, suitable for use in larger analysis pipelines.

## Usage

### Requirements

- [Cromwell](https://cromwell.readthedocs.io/), [MiniWDL](https://github.com/chanzuckerberg/miniwdl), [Sprocket](https://sprocket.bio/), or another WDL-compatible workflow executor
- Docker/Apptainer (the workflow uses `getwilds/bwa:0.7.17` container)

### Importing into Your Workflow

TBD

### Integration Examples

TBD

## Testing the Module

TBD

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:

- **Memory**: TBD
- **CPUs**: TBD

### Advanced Parameters

- Resource allocation can be tuned for different compute environments

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Sufficient memory for alignment

## Features

- **BWA-MEM alignment**: Align sequences 70bp to 1Mbp in length
- **Validation**: Built-in output validation and reporting
- **Scalable**: Configurable resource allocation
- **Robust**: Extensive error handling and cleanup
- **Compatible**: Works with multiple WDL executors

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real RNA-seq data (chromosome 22 subset for efficiency)
- Comprehensive validation of all outputs

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
