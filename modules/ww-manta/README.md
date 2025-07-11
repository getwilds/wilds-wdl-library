# ww-manta
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for structural variant calling using Manta.

## Overview

This module provides reusable WDL tasks for high-quality structural variant detection using Manta. It supports both DNA-seq and RNA-seq data and includes comprehensive validation of outputs. Manta is optimized for detecting large insertions, deletions, duplications, inversions, and translocations from sequencing data.

The module uses Manta's algorithms for precise structural variant detection and can run completely standalone with automatic test data download, or integrate with existing BAM files.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `manta_call`, `validate_outputs`
- **Workflow**: `manta_example` (demonstration workflow with automatic test data support)
- **Container**: `getwilds/manta:1.6.0`
- **Dependencies**: Integrates with `ww-testdata` module for complete workflows
- **Test Data**: Automatically downloads reference genome and BAM data when not provided

## Tasks

### `manta_call`
Calls structural variants using Manta for a single sample.

**Inputs:**
- `aligned_bam` (File): Input aligned BAM file containing reads for variant calling
- `aligned_bam_index` (File): Index file for the aligned BAM above
- `sample_name` (String): Name of the sample provided for output files
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `call_regions_bed` (File?): Optional BED file to restrict calling to specific regions
- `call_regions_index` (File?): Optional index file for the BED file
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
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-manta/ww-manta.wdl" as manta_tasks

struct MantaSample {
    String name
    File bam
    File bai
}

workflow my_sv_pipeline {
  input {
    Array[MantaSample] samples
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

### Advanced Usage Examples

**RNA-seq mode:**
```wdl
call manta_tasks.manta_call {
  input:
    is_rna = true,  # Enable RNA-specific algorithms
}
```

**Targeted calling:**
```wdl
call manta_tasks.manta_call {
  input:
    call_regions_bed = target_regions.bed,
    call_regions_index = target_regions.bed.tbi,
    # Focus calling on specific genomic regions
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-bwa**: Use aligned BAM files from BWA as input
- **ww-star**: Use RNA-seq alignments for RNA structural variant calling
- **ww-testdata**: Automatic provisioning of reference data and test samples
- **Comprehensive genomics pipelines**: Combine with SNV/indel callers for complete variant analysis

## Testing the Module

The module includes a demonstration workflow that can be tested independently:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-manta.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-manta.wdl -i inputs.json

# Using Sprocket
sprocket run ww-manta.wdl inputs.json
```

### Automatic Demo Mode

When no samples or reference files are provided, the workflow automatically:
1. Downloads reference genome data using `ww-testdata`
2. Downloads demonstration BAM data using `ww-testdata`
3. Calls structural variants using Manta
4. Validates all outputs

### Test Input Format

**Minimal input (uses automatic demo data):**
```json
{
  "manta_example.is_rna": false,
  "manta_example.cpus": 2,
  "manta_example.memory_gb": 8
}
```

**Full input (provide your own data):**
```json
{
  "manta_example.samples": [
    {
      "name": "sample1",
      "bam": "/path/to/sample1.bam",
      "bai": "/path/to/sample1.bam.bai"
    }
  ],
  "manta_example.ref_fasta": "/path/to/reference.fasta",
  "manta_example.ref_fasta_index": "/path/to/reference.fasta.fai",
  "manta_example.call_regions_bed": "/path/to/target_regions.bed",
  "manta_example.call_regions_index": "/path/to/target_regions.bed.tbi",
  "manta_example.is_rna": false,
  "manta_example.cpus": 8,
  "manta_example.memory_gb": 16
}
```

**Note**: You can mix and match - provide some inputs and let others use test data.

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:

- **Memory**: 8-16 GB recommended for most analyses, may need more for very high-coverage samples
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

### Demo Configuration

- `is_rna`: Use false for faster demo runs with DNA-seq algorithms
- Resource parameters apply to both demo and user-provided data modes

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Input BAM files must be coordinate-sorted and indexed (when providing your own data)
- Reference genome FASTA file with associated index (when providing your own data)

## Features

- **Standalone execution**: Complete workflow with automatic test data download
- **Flexible input**: Use your own data or automatic demo data
- **Comprehensive SV detection**: Detects deletions, insertions, duplications, inversions, and translocations
- **Multi-sample support**: Process multiple samples in parallel
- **RNA-seq compatibility**: Specialized mode for RNA structural variant detection
- **Validation**: Built-in output validation and reporting
- **Module integration**: Seamlessly combines with ww-testdata
- **Scalable**: Configurable resource allocation
- **Robust**: Extensive error handling and cleanup
- **Compatible**: Works with multiple WDL executors

## Output Description

- **VCF files**: Contain structural variant calls in compressed VCF format with proper INFO and FORMAT fields
- **VCF indices**: Enable rapid random access to variant regions (.tbi format)
- **Validation report**: Comprehensive validation with file integrity checks, format validation, and variant counts

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data (demonstration BAM for integration testing)
- Comprehensive validation of all outputs including VCF format validation
- Integration testing with ww-testdata modules

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to Manta usage or configuration, please refer to the documentation present in the [Manta GitHub repository](https://github.com/Illumina/manta). Please make sure to cite their work if you use Manta in your analyses.

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
