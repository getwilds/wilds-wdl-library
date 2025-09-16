# ww-bcftools
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for variant calling and manipulation using bcftools.

## Overview

This module provides reusable WDL tasks for variant calling using bcftools, a popular suite of utilities for variant calling and manipulating VCF/BCF files. The module implements bcftools' mpileup/call pipeline for variant detection from BAM files with comprehensive parameter customization.

The module is designed to be a foundational component within the WILDS ecosystem, suitable for use in larger variant calling pipelines, and includes automatic test data downloading when no inputs are provided.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `mpileup_call`, `validate_outputs`  
- **Workflow**: `bcftools_example` (demonstration workflow that executes core tasks)
- **Container**: `getwilds/bcftools:1.19`
- **Test Data**: Automatically downloads test data when no samples provided using `ww-testdata` module

## Tasks

### `mpileup_call`
Calls variants using bcftools mpileup and call pipeline with comprehensive parameter support.

**Inputs:**
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Reference genome FASTA index (.fai)
- `bam_file` (File): Input BAM file for variant calling
- `bam_index` (File): BAM index file (.bai)
- `regions_bed` (File?): Optional BED file for targeted regions
- `annotate_format` (String): FORMAT annotations (default: "FORMAT/AD,FORMAT/DP")
- `ignore_rg` (Boolean): Ignore read groups during analysis (default: true)
- `disable_baq` (Boolean): Disable BAQ computation (default: true)
- `max_depth` (Int): Maximum read depth for mpileup (default: 10000)
- `max_idepth` (Int): Maximum per-sample indel depth (default: 10000)
- `memory_gb` (Int): Memory allocation (default: 8)
- `cpu_cores` (Int): CPU cores (default: 2)

**Outputs:**
- `mpileup_vcf` (File): Compressed VCF file with called variants
- `mpileup_vcf_index` (File): Index for the VCF file (.csi)

### `validate_outputs`
Validates variant calling outputs and generates comprehensive statistics including variant counts and SNP/indel breakdowns.

**Inputs:**
- `vcf_files` (Array[File]): Array of VCF files to validate

**Outputs:**
- `report` (File): Validation summary with detailed variant statistics

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bcftools/ww-bcftools.wdl" as bcftools_tasks

struct BcftoolsSample {
    String name
    File bam
    File bai
}

workflow my_variant_calling_pipeline {
  input {
    Array[BcftoolsSample] samples
    File reference_fasta
    File reference_fasta_index
    File? regions_bed
  }
  
  scatter (sample in samples) {
    call bcftools_tasks.mpileup_call {
      input:
        bam_file = sample.bam,
        bam_index = sample.bai,
        reference_fasta = reference_fasta,
        reference_fasta_index = reference_fasta_index,
        regions_bed = regions_bed,
        max_depth = 8000,
        max_idepth = 8000
    }
  }
  
  call bcftools_tasks.validate_outputs {
    input:
      vcf_files = mpileup_call.mpileup_vcf
  }
  
  output {
    Array[File] variant_vcfs = mpileup_call.mpileup_vcf
    Array[File] vcf_indices = mpileup_call.mpileup_vcf_index
    File validation_report = validate_outputs.report
  }
}
```

### Integration Examples

This module integrates well with other WILDS components:
- **ww-testdata**: Provides reference genomes and test BAM files
- **ww-star**: Use aligned BAM files from STAR for variant calling
- **ww-bwa**: Use aligned BAM files from BWA for variant calling
- **ww-annotsv**: Annotate structural variants found in VCF outputs
- **Custom workflows**: Combine with GATK or other variant callers for consensus calling

## Testing the Module

The module includes a demonstration workflow with comprehensive testing:

The demonstration workflow automatically downloads test data and runs without requiring input files:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-bcftools.wdl

# Using miniWDL
miniwdl run ww-bcftools.wdl

# Using Sprocket
sprocket run ww-bcftools.wdl
```

The demonstration workflow (`bcftools_example`) automatically downloads test data using the `ww-testdata` module and runs variant calling without requiring any input parameters.

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:
- **Memory**: 4-8GB recommended for most applications, scale with BAM size
- **CPUs**: 1-2 cores typically sufficient for bcftools operations
- **Storage**: Sufficient space for input BAMs and output VCFs

### Variant Calling Parameters

**Depth Settings:**
- `max_depth`: Adjust based on expected coverage depth (default: 10000)
- `max_idepth`: Set higher for high-coverage applications (default: 10000)

**Quality Settings:**
- `annotate_format`: Customize FORMAT fields in output VCF (default: "FORMAT/AD,FORMAT/DP")
- `ignore_rg`: Set to false if read group information is important
- `disable_baq`: Set to false for higher sensitivity (may increase false positives)

**Targeted Analysis:**
- `regions_bed`: Use for targeted sequencing applications to focus on specific regions

## Advanced Usage

The module supports various bcftools mpileup/call parameters for fine-tuning variant calling:

```wdl
call bcftools_tasks.mpileup_call {
  input:
    bam_file = sample.bam,
    bam_index = sample.bai,
    regions_bed = my_targets_bed,
    reference_fasta = ref_fasta,
    reference_fasta_index = ref_fasta_index,
    max_depth = 1000,
    max_idepth = 1000,
    annotate_format = "FORMAT/AD,FORMAT/DP,FORMAT/SP",
    ignore_rg = false,
    disable_baq = false
}
```

## Output Format

The module generates standard VCF files with the following characteristics:
- **Compression**: bgzip-compressed (.vcf.gz)
- **Indexing**: CSI index files (.csi) for rapid access
- **FORMAT fields**: Configurable annotations (default: AD, DP)
- **Variant types**: SNPs, indels, and complex variants

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support for containerized execution
- Docker/Apptainer support
- Sufficient computational resources (8GB RAM recommended)

## Features

- **Flexible variant calling**: Configurable mpileup/call parameters
- **Targeted regions**: Support for BED file region specification  
- **Automatic test data**: Downloads test data when no inputs provided
- **Comprehensive validation**: Built-in output validation with variant statistics
- **Scalable**: Configurable resource allocation per sample
- **Compatible**: Works with multiple WDL executors
- **Integration ready**: Designed for use with other WILDS modules

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data (chromosome 1 subset for efficiency)
- Comprehensive validation of all outputs
- Cross-platform compatibility testing

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
