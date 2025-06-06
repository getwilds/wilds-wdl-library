# ww-bcftools
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for variant calling and manipulation using bcftools.

## Overview

This module provides reusable WDL tasks for variant calling and manipulation using bcftools, a popular suite of utilities for variant calling and manipulating VCF/BCF files. The module implements bcftools' mpileup/call pipeline for variant detection and includes additional tasks for consensus sequence generation and variant filtering.

The module is designed to be a foundational component within the WILDS ecosystem, suitable for use in larger variant calling pipelines.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `mpileup_call`, `validate_outputs`
- **Workflow**: `bcftools_example` (demonstration workflow that executes core tasks)
- **Container**: `getwilds/bcftools:1.19`

## Tasks

### `mpileup_call`
Calls variants using bcftools mpileup and call pipeline.

**Inputs:**
- `sample_data` (SampleInfo): Sample information including BAM files
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Reference genome FASTA index
- `regions_bed` (File?): Optional BED file for targeted regions
- `max_depth` (Int): Maximum read depth (default: 10000)
- `max_idepth` (Int): Maximum indel depth (default: 10000)
- Various additional parameters for fine-tuning

**Outputs:**
- `output_vcf` (File): Compressed VCF file with called variants
- `output_vcf_index` (File): Index for the VCF file

### `validate_outputs`
Validates variant calling outputs and generates comprehensive statistics.

**Inputs:**
- Arrays of output files from variant calling tasks

**Outputs:**
- `report` (File): Validation summary with variant statistics

## Usage as a Module

### Importing into Your Workflow

```wdl
import "path/to/ww-bcftools.wdl" as bcftools_tasks

struct SampleInfo {
    String name
    File bam
    File bai
}

workflow my_variant_calling_pipeline {
  input {
    Array[SampleInfo] samples
    File reference_fasta
    File reference_fasta_index
    File? regions_bed
  }
  
  scatter (sample in samples) {
    call bcftools_tasks.mpileup_call {
      input:
        sample_data = sample,
        reference_fasta = reference_fasta,
        reference_fasta_index = reference_fasta_index,
        regions_bed = regions_bed
    }
  }
  
  output {
    Array[File] variant_vcfs = mpileup_call.output_vcf
  }
}
```

### Integration Examples

This module integrates well with other WILDS components:
- **ww-star**: Use aligned BAM files from STAR for variant calling
- **ww-sra**: Download data, align, then call variants
- **Custom workflows**: Combine with GATK or other variant callers for consensus calling

## Testing the Module

The module includes a demonstration workflow with comprehensive testing:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-bcftools.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-bcftools.wdl -i inputs.json

# Using Sprocket
sprocket run ww-bcftools.wdl inputs.json
```

### Test Input Format

```json
{
  "bcftools_example.samples": [
    {
      "name": "sample1",
      "bam": "/path/to/sample1.bam",
      "bai": "/path/to/sample1.bam.bai"
    }
  ],
  "bcftools_example.reference_genome": {
    "name": "hg38",
    "fasta": "/path/to/genome.fasta",
    "fasta_index": "/path/to/genome.fasta.fai"
  },
  "bcftools_example.regions_bed": "/path/to/regions.bed",
  "bcftools_example.max_depth": 10000,
  "bcftools_example.memory_gb": 8,
  "bcftools_example.cpu_cores": 2
}
```

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:
- **Memory**: 4-8GB recommended for most applications
- **CPUs**: 1-2 cores typically sufficient for bcftools operations
- **Storage**: Sufficient space for input BAMs and output VCFs

### Variant Calling Parameters

- `max_depth`: Adjust based on expected coverage depth
- `max_idepth`: Set higher for high-coverage applications
- `annotate_format`: Customize FORMAT fields in output VCF
- `regions_bed`: Use for targeted sequencing applications

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Input BAM files must be sorted and indexed
- Reference genome with FASTA index

## Features

- **Variant calling**: Configurable mpileup/call parameters
- **Targeted regions**: Support for BED file region specification
- **Validation**: Built-in output validation and statistics
- **Scalable**: Configurable resource allocation
- **Compatible**: Works with multiple WDL executors

## Advanced Usage

The module supports various bcftools mpileup/call parameters for fine-tuning variant calling:

```wdl
call bcftools_tasks.mpileup_call {
  input:
    sample_data = my_sample,
    reference_fasta = ref_fasta,
    max_depth = 1000,
    max_idepth = 1000,
    annotate_format = "FORMAT/AD,FORMAT/DP,FORMAT/SP",
    ignore_rg = false,
    disable_baq = false
}
```

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data (chromosome 22 subset for efficiency)
- Comprehensive validation of all outputs

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
