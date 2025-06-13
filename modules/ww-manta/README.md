# ww-manta
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for structural variant calling using Manta.

## Overview

This module provides reusable WDL tasks for high-quality structural variant detection using Manta. It supports both DNA-seq and RNA-seq data and includes comprehensive validation of outputs. Manta is optimized for detecting large insertions, deletions, duplications, inversions, and translocations from sequencing data.

The module is designed to be a foundational component within the WILDS ecosystem, suitable for use in larger genomics analysis pipelines.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `manta_call`, `validate_outputs`
- **Workflow**: `manta_example` (demonstration workflow that executes all tasks)
- **Container**: `getwilds/manta:1.6.0`

## Tasks

### `manta_call`
Calls structural variants using Manta for a single sample.

**Inputs:**
- `aligned_bam` (File): Input aligned BAM file containing reads for variant calling
- `aligned_bam_index` (File): Index file for the aligned BAM above
- `sample_name` (String): Name of the sample provided for output files
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `call_regions_bed` (File, optional): BED file to restrict calling to specific regions
- `call_regions_index` (File, optional): Index file for the optional BED file
- `is_rna` (Boolean): Flag for RNA-seq mode (default: false)
- `cpu_cores` (Int): Number of CPU cores (default: 8)
- `memory_gb` (Int): Memory allocation in GB (default: 16)

**Outputs:**
- `vcf` (File): Compressed VCF file with structural variant calls
- `vcf_index` (File): Index file for the VCF

### `validate_outputs`
Validates Manta outputs and generates a comprehensive report.

**Inputs:**
- `vcf_files` (Array[File]): Array of VCF files to validate
- `vcf_index_files` (Array[File]): Array of VCF index files to validate

**Outputs:**
- `report` (File): Validation summary with structural variant statistics

## Usage as a Module

### Importing into Your Workflow

```wdl
import "path/to/ww-manta.wdl" as manta_tasks

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
    call manta_tasks.manta_call {
      input:
        aligned_bam = sample.bam,
        aligned_bam_index = sample.bai,
        sample_name = sample.name,
        reference_fasta = reference_fasta,
        reference_fasta_index = reference_fasta_index
    }
  }
  
  call manta_tasks.validate_outputs {
    input:
      vcf_files = manta_call.vcf,
      vcf_index_files = manta_call.vcf_index
  }
  
  output {
    Array[File] sv_vcfs = manta_call.vcf
    Array[File] sv_vcf_indexes = manta_call.vcf_index
    File validation_report = validate_outputs.report
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-bwa**: Use aligned BAM files from BWA as input
- **ww-star**: Use RNA-seq alignments for RNA structural variant calling
- **Comprehensive genomics pipelines**: Combine with SNV/indel callers for complete variant analysis

## Testing the Module

The module includes a demonstration workflow with comprehensive testing:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-manta.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-manta.wdl -i inputs.json

# Using Sprocket
sprocket run ww-manta.wdl inputs.json
```

### Test Input Format

```json
{
  "manta_example.samples": [
    {
      "name": "sample1",
      "bam": "/path/to/sample1.bam",
      "bai": "/path/to/sample1.bam.bai"
    }
  ],
  "manta_example.reference_genome": {
    "name": "hg38",
    "fasta": "/path/to/genome.fasta",
    "fasta_index": "/path/to/genome.fasta.fai"
  },
  "manta_example.cpus": 8,
  "manta_example.memory_gb": 16
}
```

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:

- **Memory**: 16GB recommended for most analyses, may need more for very high-coverage samples
- **CPUs**: 8 cores recommended, can be adjusted based on available resources
- **Storage**: Ensure sufficient space for intermediate files (typically 2-3x input BAM size)

### Data Type Considerations

- **DNA-seq**: Use default settings (`is_rna = false`)
- **RNA-seq**: Set `is_rna = true` to enable RNA-specific algorithms that account for splicing
- **Targeted panels**: Use `call_regions_bed` to restrict analysis to regions of interest

### Advanced Parameters

- **call_regions_bed**: Restrict calling to specific genomic regions for targeted sequencing
- **is_rna**: Enables RNA-seq mode with splice-aware algorithms
- Resource allocation can be tuned for different compute environments

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Input BAM files must be coordinate-sorted and indexed
- Reference genome FASTA file with associated index

## Features

- **Comprehensive SV detection**: Detects deletions, insertions, duplications, inversions, and translocations
- **Multi-sample support**: Process multiple samples in parallel
- **RNA-seq compatibility**: Specialized mode for RNA structural variant detection
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
