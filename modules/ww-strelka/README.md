# ww-strelka

A WILDS WDL module for germline and somatic variant calling using Strelka, a fast and accurate small variant caller optimized for analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs.

## Overview

Strelka is a small variant caller designed to detect single nucleotide variants (SNVs) and small insertions/deletions (indels) from Illumina sequencing data. This module provides both germline and somatic variant calling capabilities within a standardized WILDS framework.

### Key Features

- **Dual calling modes**: Supports both germline and somatic variant calling
- **Zero-configuration testing**: Test workflow requires no input parameters
- **Exome optimization**: Specialized settings for exome sequencing analysis
- **Parallel processing**: Multi-sample support with configurable resource allocation
- **Comprehensive validation**: Built-in output validation and quality reporting
- **Flexible targeting**: Tasks support optional BED file targeting for custom workflows
- **Module integration**: Seamlessly integrates with other WILDS modules

## Test Workflow

### No Input Required

The `strelka_example` test workflow requires no input parameters and automatically:

1. Downloads test reference genome data via `ww-testdata`
2. Downloads test tumor BAM data via `ww-testdata`
3. Downloads test normal BAM data via `ww-testdata`
4. Performs both germline and somatic variant calling
5. Validates all outputs and generates comprehensive reports

### Hardcoded Test Settings

- **Mode**: Both germline and somatic calling
- **Exome**: Disabled (whole genome mode)
- **Resources**: 4 CPUs, 8GB RAM per task
- **Samples**: Single tumor/normal pair from test data
- **Target regions**: None (genome-wide calling)

### StrelkaSample Structure

```wdl
struct StrelkaSample {
    String name      # Sample identifier
    File bam         # Coordinate-sorted BAM file
    File bai         # BAM index file
}
```

## Outputs

| Output | Type | Description |
|--------|------|-------------|
| `germline_vcfs` | Array[File] | Germline variant calls in compressed VCF format |
| `germline_vcf_indices` | Array[File] | Index files for germline VCFs |
| `somatic_snvs_vcfs` | Array[File] | Somatic SNV calls in compressed VCF format |
| `somatic_indels_vcfs` | Array[File] | Somatic indel calls in compressed VCF format |
| `somatic_snvs_vcf_indices` | Array[File] | Index files for somatic SNV VCFs |
| `somatic_indels_vcf_indices` | Array[File] | Index files for somatic indel VCFs |
| `validation_report` | File | Combined validation report with file checks and statistics |

## Usage Examples

### Germline Variant Calling

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-strelka/ww-strelka.wdl" as strelka

workflow my_germline_analysis {
  input {
    Array[File] bam_files
    Array[File] bai_files
    File reference_fasta
    File reference_index
  }
  
  # Create sample array with zip (bam files on the 'left', their indices on the 'right')
  Array[Pairs] sample_pairs = zip(bam_files, bai_files)
  scatter (pair in sample_pairs) {
    String sample_name = basename(pair.left, ".bam")
    call strelka.strelka_germline { input:
      sample_name = sample_name,
      bam = pair.left,
      bai = pair.right,
      ref_fasta = reference_fasta,
      ref_fasta_index = reference_index,
      is_exome = true  # Enable exome optimizations
    }
  }
}
```

### Somatic Variant Calling

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-strelka/ww-strelka.wdl" as strelka

workflow my_somatic_analysis {
  input {
    File tumor_bam
    File tumor_bai
    File normal_bam
    File normal_bai
    File reference_fasta
    File reference_index
  }
  
  call strelka.strelka_somatic { input:
    tumor_sample_name = "tumor",
    tumor_bam = tumor_bam,
    tumor_bai = tumor_bai,
    normal_sample_name = "normal",
    normal_bam = normal_bam,
    normal_bai = normal_bai,
    ref_fasta = reference_fasta,
    ref_fasta_index = reference_index
  }
}
```

### Running the Test Workflow

```bash
# No input file needed for test workflow
miniwdl run ww-strelka.wdl

# Or with other executors
java -jar cromwell.jar run ww-strelka.wdl
sprocket run ww-strelka.wdl
```

## Workflow Tasks

### strelka_germline

Performs germline variant calling on individual samples.

**Inputs:**
- Sample name
- Sample BAM file and index
- Reference genome FASTA and index
- Optional target regions BED file
- Exome sequencing flag

**Outputs:**
- Compressed germline variants VCF
- VCF index file

### strelka_somatic

Performs somatic variant calling on tumor/normal pairs.

**Inputs:**
- Tumor and normal sample names
- Tumor and normal BAM files with indices
- Reference genome FASTA and index
- Optional target regions BED file
- Exome sequencing flag

**Outputs:**
- Compressed somatic SNVs VCF
- Compressed somatic indels VCF
- VCF index files for both outputs

### Validation Tasks

The module includes comprehensive validation tasks that:
- Verify all output files are present and non-empty
- Validate VCF format compliance
- Count variants and generate summary statistics
- Combine reports from germline and somatic analyses

## Requirements

- **WDL executor**: Cromwell, miniWDL, Sprocket, or compatible executor
- **Docker/Apptainer**: Container runtime for reproducible execution
- **Input data**: Coordinate-sorted BAM files with indices (when providing your own data)
- **Reference genome**: FASTA file with index (when providing your own data)
- **Resources**: Minimum 4 CPU cores and 8GB RAM (configurable)

## Performance Considerations

- **CPU scaling**: Strelka scales effectively with increased CPU allocation
- **Memory requirements**: 8GB is typically sufficient for most analyses
- **Runtime**: Whole genome germline calling typically completes in 1-3 hours
- **Somatic calling**: Generally faster than germline due to targeted approach
- **Exome data**: Significantly faster processing with `is_exome=true`

## Output Description

### Germline Calling
- **variants.vcf.gz**: All germline variant calls with quality scores and filters
- **variants.vcf.gz.tbi**: Tabix index for rapid region-based queries

### Somatic Calling
- **somatic.snvs.vcf.gz**: Somatic single nucleotide variants
- **somatic.indels.vcf.gz**: Somatic insertions and deletions  
- **\*.tbi files**: Tabix indices for all VCF outputs

### Validation Reports
- File integrity verification
- Variant count summaries
- Format validation results
- Overall workflow status

## Algorithm Details

Strelka uses probabilistic models to:
- Distinguish true variants from sequencing artifacts
- Optimize sensitivity and specificity for different variant types
- Handle low-frequency somatic variants in tumor samples
- Account for normal tissue contamination in tumor samples

The tool is particularly effective for:
- High-confidence germline variant discovery
- Somatic mutation detection in tumor/normal pairs
- Analysis of targeted sequencing panels
- Whole genome and exome sequencing data

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data for integration testing
- Comprehensive validation of all outputs
- Integration testing with ww-testdata modules

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to Strelka usage, configuration, or algorithm details, please refer to the [Strelka GitHub repository](https://github.com/Illumina/strelka) and associated documentation. Please make sure to cite Strelka in your publications if you use this module.

## Citation

If you use this module in your research, please cite:

**Strelka:**
Kim, S., Scheffler, K., Halpern, A.L. et al. Strelka2: fast and accurate calling of germline and somatic variants. Nat Methods 15, 591–594 (2018). https://doi.org/10.1038/s41592-018-0051-x

**WILDS:**
Please also acknowledge the WILDS project in your publications.

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
