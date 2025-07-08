# ww-bcftools
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for variant calling and manipulation using bcftools.

## Overview

This module provides reusable WDL tasks for variant calling and manipulation using bcftools, a popular suite of utilities for variant calling and manipulating VCF/BCF files. The module implements bcftools' mpileup/call pipeline for variant detection and includes comprehensive validation tasks.

The module is designed to be a foundational component within the WILDS ecosystem, suitable for use in larger variant calling pipelines. It can run completely standalone with automatic test data download and alignment, or integrate with existing BAM files.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `mpileup_call`, `validate_outputs`
- **Workflow**: `bcftools_example` (demonstration workflow with automatic test data support)
- **Container**: `getwilds/bcftools:1.19`
- **Dependencies**: Integrates with `ww-sra`, `ww-bwa`, and `ww-testdata` modules for complete workflows
- **Test Data**: Automatically downloads reference genome and SRA data when not provided

## Tasks

### `mpileup_call`
Calls variants using bcftools mpileup and call pipeline.

**Inputs:**
- `bam_file` (File): Input BAM file for the sample
- `bam_index` (File): Index file for the input BAM
- `reference_fasta` (File): Reference genome FASTA file  
- `reference_fasta_index` (File): Reference genome FASTA index
- `regions_bed` (File?): Optional BED file for targeted regions
- `annotate_format` (String): FORMAT annotations to add (default: "FORMAT/AD,FORMAT/DP")
- `ignore_rg` (Boolean): Ignore read groups during analysis (default: true)
- `disable_baq` (Boolean): Disable BAQ computation (default: true)
- `max_depth` (Int): Maximum read depth (default: 10000)
- `max_idepth` (Int): Maximum indel depth (default: 10000)
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores (default: 2)

**Outputs:**
- `mpileup_vcf` (File): Compressed VCF file with called variants
- `mpileup_vcf_index` (File): Index for the VCF file

### `validate_outputs`
Validates variant calling outputs and generates comprehensive statistics.

**Inputs:**
- `vcf_files` (Array[File]): Array of VCF files to validate

**Outputs:**
- `report` (File): Validation summary with variant statistics

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
        regions_bed = regions_bed
    }
  }
  
  output {
    Array[File] variant_vcfs = mpileup_call.mpileup_vcf
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-star**: Use aligned BAM files from STAR for variant calling
- **ww-bwa**: Use aligned BAM files from BWA for variant calling
- **ww-sra**: Download data, align, then call variants (built into demo workflow)
- **Custom workflows**: Combine with GATK or other variant callers for consensus calling

## Testing the Module

The module includes a demonstration workflow with comprehensive testing and **automatic test data download**:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-bcftools.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-bcftools.wdl -i inputs.json

# Using Sprocket
sprocket run ww-bcftools.wdl inputs.json
```

### Automatic Demo Mode

When no samples or reference files are provided, the workflow automatically:
1. Downloads reference genome data using `ww-testdata`
2. Downloads SRA data (default: ERR1258306) using `ww-sra`  
3. Builds BWA index using `ww-bwa`
4. Aligns reads using `ww-bwa`
5. Calls variants using `bcftools`

### Test Input Format

**Minimal input (uses automatic demo data):**
```json
{
  "bcftools_example.demo_sra_id": "ERR1258306",
  "bcftools_example.max_depth": 10000,
  "bcftools_example.memory_gb": 8,
  "bcftools_example.cpu_cores": 2
}
```

**Full input (provide your own data):**
```json
{
  "bcftools_example.samples": [
    {
      "name": "sample1", 
      "bam": "/path/to/sample1.bam",
      "bai": "/path/to/sample1.bam.bai"
    }
  ],
  "bcftools_example.ref_fasta": "/path/to/genome.fasta",
  "bcftools_example.ref_fasta_index": "/path/to/genome.fasta.fai",
  "bcftools_example.regions_bed": "/path/to/regions.bed",
  "bcftools_example.max_depth": 10000,
  "bcftools_example.max_idepth": 10000,
  "bcftools_example.memory_gb": 8,
  "bcftools_example.cpu_cores": 2
}
```

**Note**: You can mix and match - provide some inputs and let others use test data.

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
- `ignore_rg`: Set to false if read group information is important
- `disable_baq`: Set to false for more sensitive calling
- `regions_bed`: Use for targeted sequencing applications

### Demo Configuration

- `demo_sra_id`: Change to use different SRA sample for testing
- Resource parameters apply to both demo and user-provided data modes

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Input BAM files must be sorted and indexed (when providing your own data)
- Reference genome with FASTA index (when providing your own data)

## Features

- **Standalone execution**: Complete workflow with automatic test data download
- **Flexible input**: Use your own data or automatic demo data
- **Variant calling**: Configurable mpileup/call parameters
- **Targeted regions**: Support for BED file region specification  
- **Validation**: Built-in output validation and statistics
- **Module integration**: Seamlessly combines with ww-sra, ww-bwa, and ww-testdata
- **Scalable**: Configurable resource allocation
- **Compatible**: Works with multiple WDL executors

## Advanced Usage

The module supports various bcftools mpileup/call parameters for fine-tuning variant calling:

```wdl
call bcftools_tasks.mpileup_call {
  input:
    bam_file = my_sample.bam,
    bam_index = my_sample.bai,
    reference_fasta = ref_fasta,
    reference_fasta_index = ref_fasta_index,
    regions_bed = my_bedfile,
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
- Real sequencing data (SRA sample ERR1258306 for integration testing)
- Comprehensive validation of all outputs
- Integration testing with ww-sra, ww-bwa, and ww-testdata modules

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
