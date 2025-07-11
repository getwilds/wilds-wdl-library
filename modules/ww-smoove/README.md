# ww-smoove
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for structural variant calling using Smoove.

## Overview

This module provides reusable WDL tasks for structural variant (SV) calling using Smoove, a streamlined wrapper around lumpy-sv that simplifies structural variant detection. Smoove excels at calling deletions, duplications, inversions, and translocations from whole-genome sequencing data with high sensitivity and specificity.

The module can run completely standalone with automatic test data download, or integrate with existing BAM files for production analyses.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `smoove_call`, `validate_outputs`
- **Workflow**: `smoove_example` (demonstration workflow with automatic test data support)
- **Container**: `brentp/smoove:latest`
- **Dependencies**: Integrates with `ww-testdata` module for complete workflows
- **Test Data**: Automatically downloads reference genome and BAM data when not provided

## Tasks

### `smoove_call`
Calls structural variants using Smoove for a single sample.

**Inputs:**
- `aligned_bam` (File): Input aligned BAM file containing reads for variant calling
- `aligned_bam_index` (File): Index file for the aligned BAM above
- `sample_name` (String): Name of the sample provided for output files
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `target_regions_bed` (File?): Optional BED file defining target regions to include in final output
- `exclude_bed` (File?): Optional BED file defining regions to exclude from calling
- `exclude_chroms` (String?): Optional comma-separated list of chromosomes to exclude
- `cpu_cores` (Int): Number of CPU cores to use (default: 8)
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

struct SmooveSample {
    String name
    File bam
    File bai
}

workflow my_sv_pipeline {
  input {
    Array[SmooveSample] samples
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

### Advanced Usage Examples

**Excluding problematic regions:**
```wdl
call smoove_tasks.smoove_call {
  input:
    exclude_bed = centromeres_telomeres.bed,
    # Recommended for whole-genome analysis
}
```

**Targeting specific regions:**
```wdl
call smoove_tasks.smoove_call {
  input:
    target_regions_bed = exons.bed,
    # Focus output on specific genomic regions
}
```

**Excluding specific chromosomes:**
```wdl
call smoove_tasks.smoove_call {
  input:
    exclude_chroms = "chrY,chrM",
    # Exclude Y chromosome and mitochondrial DNA
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-bwa**: Use aligned BAM files from BWA as input
- **ww-star**: Use RNA-seq alignments (though Smoove is optimized for DNA-seq)
- **ww-testdata**: Automatic provisioning of reference data and test samples
- **Custom workflows**: Foundation for any structural variant analysis pipeline

## Testing the Module

The module includes a demonstration workflow that can be tested independently:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-smoove.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-smoove.wdl -i inputs.json

# Using Sprocket
sprocket run ww-smoove.wdl inputs.json
```

### Automatic Demo Mode

When no samples or reference files are provided, the workflow automatically:
1. Downloads reference genome data using `ww-testdata`
2. Downloads demonstration BAM data using `ww-testdata`
3. Calls structural variants using Smoove
4. Validates all outputs

### Test Input Format

**Minimal input (uses automatic demo data):**
```json
{
  "smoove_example.cpus": 2,
  "smoove_example.memory_gb": 8
}
```

**Full input (provide your own data):**
```json
{
  "smoove_example.samples": [
    {
      "name": "sample1",
      "bam": "/path/to/sample1.bam",
      "bai": "/path/to/sample1.bam.bai"
    }
  ],
  "smoove_example.ref_fasta": "/path/to/reference.fasta",
  "smoove_example.ref_fasta_index": "/path/to/reference.fasta.fai",
  "smoove_example.exclude_bed": "/path/to/exclude_regions.bed",
  "smoove_example.include_bed": "/path/to/target_regions.bed",
  "smoove_example.cpus": 8,
  "smoove_example.memory_gb": 16
}
```

**Note**: You can mix and match - provide some inputs and let others use test data.

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:

- **Memory**: 8-16 GB recommended for most whole-genome analyses
- **CPUs**: 8 cores recommended; Smoove benefits from multi-threading
- **Storage**: Ensure sufficient space for intermediate files and output VCFs

### Region Filtering Options

- **exclude_bed**: Skip problematic regions (highly recommended for whole-genome analysis)
- **target_regions_bed**: Filter final output to include only variants overlapping target regions
- **exclude_chroms**: Exclude specific chromosomes from analysis
- Common exclusions: centromeres, telomeres, repetitive elements, sex chromosomes

### Advanced Parameters

- **target_regions_bed**: Applied as post-processing filter on final VCF output
- **exclude_bed**: Applied during variant calling to avoid problematic regions
- **exclude_chroms**: Useful for excluding sex chromosomes or mitochondrial DNA

### Demo Configuration

- Resource parameters apply to both demo and user-provided data modes

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Input BAM files must be coordinate-sorted and indexed (when providing your own data)
- Reference genome FASTA file with associated index (when providing your own data)

## Features

- **Standalone execution**: Complete workflow with automatic test data download
- **Flexible input**: Use your own data or automatic demo data
- **Streamlined SV calling**: Simplified interface to lumpy-sv algorithms
- **Multi-sample support**: Process multiple samples in parallel
- **Flexible filtering**: Multiple options for region inclusion/exclusion
- **Validation**: Built-in output validation and reporting
- **Module integration**: Seamlessly combines with ww-testdata
- **Robust**: Extensive error handling and cleanup
- **Compatible**: Works with multiple WDL executors

## Output Description

- **VCF files**: Contain structural variant calls in compressed VCF format with lumpy-sv annotations
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

For questions specific to Smoove usage or configuration, please refer to the documentation present in the [Smoove GitHub repository](https://github.com/brentp/smoove). Please make sure to cite their work if you use Smoove in your analyses.

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
