# ww-gatk
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for GATK variant calling and analysis tasks with automated parallelization.

## Overview

This module provides comprehensive variant calling capabilities using GATK (Genome Analysis Toolkit), supporting both germline and somatic variant detection workflows with **automatic interval-based parallelization** for optimal performance on whole genome sequencing (WGS) data. It includes base quality score recalibration (BQSR), germline variant calling with HaplotypeCaller, somatic variant calling with Mutect2, and quality metrics collection. The module can run completely standalone with automatic test data download or integrate with existing BAM files and reference data.

The module implements GATK best practices for variant calling, including proper base quality recalibration using known variant sites, **automatic interval splitting for parallelization**, and provides comprehensive validation and reporting for quality assurance.

## Key Features

- **Automatic parallelization**: Intelligently splits genome into optimal intervals for parallel processing
- **Significant speedup**: Near-linear performance scaling for HaplotypeCaller and Mutect2 on WGS data
- **Smart interval management**: Uses GATK SplitIntervals to create balanced computational chunks
- **Seamless merging**: Automatically combines parallel results into final VCF files
- **Zero configuration**: Works out-of-the-box with sensible defaults for parallelization

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `base_recalibrator`, `haplotype_caller`, `mutect2`, `split_intervals`, `merge_vcfs`, `merge_mutect_stats`, `create_sequence_dictionary`, `collect_wgs_metrics`, `validate_outputs`
- **Workflow**: `gatk_example` (demonstration workflow with automatic test data support and parallelization)
- **Container**: `getwilds/gatk:4.6.1.0`
- **Dependencies**: Integrates with `ww-testdata` module for automatic reference genome and variant database downloads
- **Test Data**: Automatically downloads reference genome, dbSNP, known indels, gnomAD, and aligned BAM data when not provided

## Tasks

### Core Variant Calling Tasks

### `base_recalibrator`
Performs Base Quality Score Recalibration (BQSR) to improve the accuracy of base quality scores.

**Inputs:**
- `aligned_bam` (File): Input aligned BAM file to be recalibrated
- `aligned_bam_index` (File): Index file for the input BAM
- `intervals` (File?): Optional interval list file defining target regions
- `dbsnp_vcf` (File): dbSNP VCF file for known variant sites
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `reference_dict` (File): Reference genome sequence dictionary
- `known_indels_sites_vcfs` (Array[File]): Array of VCF files with known indel sites
- `output_basename` (String): Base name for output files
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `recalibrated_bam` (File): BAM file with recalibrated base quality scores
- `recalibrated_bai` (File): Index file for the recalibrated BAM
- `recalibration_report` (File): Base recalibration report table

### `haplotype_caller`
Calls germline variants using GATK HaplotypeCaller.

**Inputs:**
- `bam` (File): Input aligned BAM file
- `bam_index` (File): Index file for the input BAM
- `intervals` (File?): Optional interval list file defining target regions
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `reference_dict` (File): Reference genome sequence dictionary
- `dbsnp_vcf` (File): dbSNP VCF file for variant annotation
- `output_basename` (String): Base name for output files
- `memory_gb` (Int): Memory allocation in GB (default: 16)
- `cpu_cores` (Int): Number of CPU cores to use (default: 4)

**Outputs:**
- `vcf` (File): Compressed VCF file containing germline variant calls
- `vcf_index` (File): Index file for the VCF output

### `mutect2`
Calls somatic variants using GATK Mutect2 in tumor-only mode with filtering.

**Inputs:**
- `bam` (File): Input aligned BAM file
- `bam_index` (File): Index file for the input BAM
- `intervals` (File?): Optional interval list file defining target regions
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `reference_dict` (File): Reference genome sequence dictionary
- `gnomad_vcf` (File): gnomAD population allele frequency VCF for germline resource
- `output_basename` (String): Base name for output files
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `vcf` (File): Compressed VCF file containing filtered somatic variant calls
- `vcf_index` (File): Index file for the Mutect2 VCF output
- `stats_file` (File): Mutect2 statistics file
- `f1r2_counts` (File): F1R2 counts for filtering

### Parallelization and Utility Tasks

### `split_intervals`
Automatically splits genome intervals into optimal chunks for parallel processing using GATK SplitIntervals.

**Inputs:**
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `reference_dict` (File): Reference genome sequence dictionary
- `intervals` (File?): Optional interval list file defining target regions to split
- `scatter_count` (Int): Number of interval files to create (default: 24)
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `interval_files` (Array[File]): Array of interval files optimized for parallel processing

### `merge_vcfs`
Merges multiple VCF files from parallel processing into a single consolidated VCF.

**Inputs:**
- `vcfs` (Array[File]): Array of VCF files to merge
- `vcf_indices` (Array[File]): Array of VCF index files
- `output_basename` (String): Base name for output files
- `reference_dict` (File): Reference sequence dictionary
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `merged_vcf` (File): Merged VCF file
- `merged_vcf_index` (File): Index for merged VCF file

### `merge_mutect_stats`
Merges Mutect2 statistics files from parallel processing.

**Inputs:**
- `stats` (Array[File]): Array of Mutect2 stats files to merge
- `output_basename` (String): Base name for output files
- `memory_gb` (Int): Memory allocation in GB (default: 4)
- `cpu_cores` (Int): Number of CPU cores to use (default: 1)

**Outputs:**
- `merged_stats` (File): Merged Mutect2 statistics file

### Supporting Tasks

### `create_sequence_dictionary`
Creates a sequence dictionary file from a reference FASTA file.

**Inputs:**
- `reference_fasta` (File): Reference genome FASTA file
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `sequence_dict` (File): Sequence dictionary file (.dict) for the reference genome

### `collect_wgs_metrics`
Collects whole genome sequencing metrics using GATK CollectWgsMetrics.

**Inputs:**
- `bam` (File): Input aligned BAM file
- `bam_index` (File): Index file for the input BAM
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `intervals` (File?): Optional interval list file defining target regions
- `output_basename` (String): Base name for output files
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)
- `minimum_mapping_quality` (Int): Minimum mapping quality for reads (default: 20)
- `minimum_base_quality` (Int): Minimum base quality for bases (default: 20)
- `coverage_cap` (Int): Maximum coverage depth to analyze (default: 250)

**Outputs:**
- `metrics_file` (File): Comprehensive WGS metrics file with coverage and quality statistics

### `validate_outputs`
Validates GATK outputs and generates comprehensive statistics report.

**Inputs:**
- `recalibrated_bams` (Array[File]): Array of recalibrated BAM files
- `recalibrated_bais` (Array[File]): Array of BAM index files
- `haplotype_vcfs` (Array[File]): Array of HaplotypeCaller VCF files
- `mutect2_vcfs` (Array[File]): Array of Mutect2 VCF files
- `wgs_metrics` (Array[File]): Array of WGS metrics files

**Outputs:**
- `report` (File): Validation summary with file checks and basic statistics

## Workflow Parallelization Configuration

**Default behavior:**
- Splits genome into 24 balanced intervals (optimal for human genome)
- Each interval runs HaplotypeCaller and Mutect2 in parallel
- Automatically merges results into final VCF files
- Near-linear speedup for WGS analysis

## Testing the Module

The module includes a demonstration workflow that runs with minimal inputs:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-gatk.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-gatk.wdl -i inputs.json

# Using Sprocket
sprocket run ww-gatk.wdl inputs.json
```

### Test Input Format

The demonstration workflow can run with an empty inputs file (`{}`) and will automatically download test data, or you can provide your own data:

```json
{
  "gatk_example.samples": [
    {
      "name": "sample1",
      "bam": "path/to/sample1.bam",
      "bai": "path/to/sample1.bam.bai"
    }
  ],
  "gatk_example.reference_fasta": "path/to/reference.fasta",
  "gatk_example.reference_fasta_index": "path/to/reference.fasta.fai",
  "gatk_example.dbsnp_vcf": "path/to/dbsnp.vcf.gz",
  "gatk_example.known_indels_sites_vcfs": ["path/to/known_indels.vcf.gz"],
  "gatk_example.gnomad_vcf": "path/to/gnomad.vcf.gz",
  "gatk_example.intervals": "path/to/intervals.list",
  "gatk_example.scatter_count": 24
}
```

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access for downloading test data (when using automatic mode)
- Sufficient computational resources (memory-intensive for whole genome analysis)
- Multiple CPU cores recommended for optimal parallelization benefits

## Features

- **Automatic interval-based parallelization**: Dramatically improves WGS performance
- **Complete GATK pipeline**: Base recalibration, germline and somatic variant calling, QC metrics
- **Intelligent interval splitting**: Uses GATK SplitIntervals for optimal load balancing
- **Seamless result merging**: Transparent combination of parallel results
- **Automatic test data**: Downloads reference genome, variant databases, and test BAM when not provided
- **Best practices implementation**: Follows GATK best practices for variant calling workflows
- **Comprehensive validation**: Built-in output validation and quality reporting
- **Flexible inputs**: Can run standalone or integrate with existing data
- **Resource optimization**: Configurable memory and CPU allocation for each task

## Performance Considerations

### Parallelization Benefits
- **HaplotypeCaller**: Up to 20x speedup on WGS with 24 intervals
- **Mutect2**: Similar dramatic performance improvements
- **Scalability**: Performance scales nearly linearly with available CPU cores
- **Memory efficiency**: Parallel tasks use memory more efficiently than single large jobs

### Resource Requirements
- **Memory usage**: 16-32GB RAM total for WGS analysis (distributed across parallel tasks)
- **CPU scaling**: **Significant benefit from multiple cores** - recommend 8-24+ CPUs for WGS
- **Storage requirements**: Ensure sufficient disk space for reference data, BAMs, and VCF outputs
- **Network access**: Initial runs require internet connectivity for downloading reference databases
- **Region targeting**: Using intervals files significantly reduces runtime and resource requirements

### Optimization Tips
- **Use appropriate scatter_count**: Balance between parallelization and overhead
- **Ensure adequate CPU allocation**: More cores = better performance for variant calling
- **Monitor memory usage**: Each parallel task needs sufficient memory

## Output Description

- **Recalibrated BAMs**: Quality-improved BAM files with recalibrated base scores
- **Germline VCFs**: HaplotypeCaller variant calls suitable for population genetics and clinical analysis (automatically merged from parallel processing)
- **Somatic VCFs**: Mutect2 tumor-only somatic variant calls with filtering applied (automatically merged from parallel processing)
- **QC metrics**: Comprehensive sequencing quality metrics including coverage statistics
- **Validation report**: Summary of pipeline execution with file verification and basic statistics

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Chromosome 1 subset data for efficient testing
- Comprehensive validation of all outputs and statistics
- Integration testing with `ww-testdata` module
- Parallelization testing with multiple interval configurations

## Integration Patterns

This module demonstrates several key patterns:
- **Conditional data download**: Automatic fallback to test data when inputs not provided
- **Resource management**: Coordinated memory allocation across compute-intensive tasks
- **Best practices workflow**: Implementation of GATK recommended variant calling pipeline
- **Comprehensive validation**: Quality assurance for complex multi-output workflows
- **Automatic parallelization**: Transparent interval-based parallel processing with result merging

## Extending the Module

This module can be extended by:
- Adding joint calling capabilities for population-scale analysis
- Integrating additional variant filtering and annotation tools
- Including structural variant calling with other GATK tools
- Adding quality control modules (e.g., FastQC, MultiQC integration)
- Customizing parallelization strategies for specific use cases

## Related WILDS Components

- **ww-testdata module**: Automatic reference data and test sample downloads
- **ww-bwa module**: Upstream alignment for generating input BAM files
- **ww-annotsv module**: Downstream structural variant annotation
- **Other variant calling modules**: Complementary variant detection approaches

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to GATK usage or configuration, please refer to the [GATK documentation](https://gatk.broadinstitute.org/hc/en-us) and [GATK forums](https://gatk.broadinstitute.org/hc/en-us/community/topics). Please make sure to cite GATK in your analyses:

Van der Auwera GA & O'Connor BD. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media.

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
