# ww-smoove
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for structural variant calling using Smoove.

## Overview

This module provides reusable WDL tasks for high-quality structural variant detection using Smoove. Smoove is a comprehensive structural variant caller that combines multiple detection algorithms (including Lumpy) to provide robust SV detection from short-read sequencing data. It is particularly effective for detecting deletions, duplications, inversions, and translocations.

The module is designed to be a foundational component within the WILDS ecosystem, suitable for use in larger genomics analysis pipelines.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `smoove_call`, `validate_outputs`
- **Workflow**: `smoove_example` (demonstration workflow that executes all tasks)
- **Container**: `brentp/smoove:latest`

## Tasks

### `smoove_call`
Calls structural variants using Smoove for a single sample.

**Inputs:**
- `aligned_bam` (File): Input aligned BAM file containing reads for variant calling
- `aligned_bam_index` (File): Index file for the aligned BAM above
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `sample_name` (String): Name of the sample provided for output files
- `target_regions_bed` (File, optional): BED file defining regions to include for calling
- `exclude_bed` (File, optional): BED file defining regions to exclude from calling
- `exclude_chromosomes` (Array[String], optional): List of chromosomes to exclude from calling
- `cpu_cores` (Int): Number of CPU cores (default: 8)
- `memory_gb` (Int): Memory allocation in GB (default: 16)

**Outputs:**
- `vcf` (File): Compressed VCF file with structural variant calls
- `vcf_index` (File): Index file for the VCF

### `validate_outputs`
Validates Smoove outputs and generates a comprehensive report.

**Inputs:**
- `vcf_files` (Array[File]): Array of VCF files to validate
- `vcf_index_files` (Array[File]): Array of VCF index files to validate

**Outputs:**
- `report` (File): Validation summary with structural variant statistics

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-smoove/ww-smoove.wdl" as smoove_tasks

struct SampleInfo {
    String name
    File bam
    File bai
}

workflow my_sv_pipeline {
  input {
    Array[SampleInfo] samples
    File reference_fasta
    File reference_fasta_index
  }
  
  scatter (sample in samples) {
    call smoove_tasks.smoove_call {
      input:
        aligned_bam = sample.bam,
        aligned_bam_index = sample.bai,
        sample_name = sample.name,
        reference_fasta = reference_fasta,
        reference_fasta_index = reference_fasta_index
    }
  }
  
  call smoove_tasks.validate_outputs {
    input:
      vcf_files = smoove_call.vcf,
      vcf_index_files = smoove_call.vcf_index
  }
  
  output {
    Array[File] sv_vcfs = smoove_call.vcf
    Array[File] sv_vcf_indexes = smoove_call.vcf_index
    File validation_report = validate_outputs.report
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-bwa**: Use aligned BAM files from BWA as input
- **ww-star**: Use RNA-seq alignments for structural variant calling
- **ww-manta**: Combine with Manta for consensus SV calling
- **Comprehensive genomics pipelines**: Combine with SNV/indel callers for complete variant analysis

## Testing the Module

The module includes a demonstration workflow with comprehensive testing:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-smoove.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-smoove.wdl -i inputs.json

# Using Sprocket
sprocket run ww-smoove.wdl inputs.json
```

### Test Input Format

```json
{
  "smoove_example.samples": [
    {
      "name": "sample1",
      "bam": "/path/to/sample1.sorted.bam",
      "bai": "/path/to/sample1.sorted.bam.bai"
    }
  ],
  "smoove_example.reference_genome": {
    "name": "hg38",
    "fasta": "/path/to/genome.fasta",
    "fasta_index": "/path/to/genome.fasta.fai"
  },
  "smoove_example.cpus": 8,
  "smoove_example.memory_gb": 16
}
```

## Configuration Guidelines

### Data Type Considerations

- **Whole genome sequencing**: Use default settings for comprehensive SV detection
- **Exome sequencing**: Consider using `target_regions_bed` to focus on exonic regions
- **Targeted panels**: Use `target_regions_bed` to restrict analysis to regions of interest

### Advanced Parameters

- **target_regions_bed**: Restrict calling to specific genomic regions for targeted sequencing
- **exclude_bed**: Exclude problematic regions (centromeres, telomeres, repetitive sequences)
- **exclude_chroms**: Exclude specific chromosome(s)
- Resource allocation can be tuned for different compute environments

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Input BAM files must be coordinate-sorted and indexed
- Reference genome FASTA file with associated index

## Features

- **Comprehensive SV detection**: Detects deletions, insertions, duplications, inversions, and translocations
- **Multi-sample support**: Process multiple samples in parallel
- **Targeted analysis**: Support for include/exclude BED files and entire chromosomes
- **Validation**: Built-in output validation and reporting
- **Scalable**: Configurable resource allocation
- **Robust**: Extensive error handling and cleanup
- **Compatible**: Works with multiple WDL executors

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data for validation
- Comprehensive validation of all outputs

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
