# ww-gatk
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for GATK variant calling and analysis tasks with automated parallelization.

## Overview

This module provides comprehensive variant calling capabilities using GATK (Genome Analysis Toolkit), supporting both germline and somatic variant detection workflows with **automatic interval-based parallelization** for optimal performance on whole genome sequencing (WGS) data. It includes duplicate marking, base quality score recalibration (BQSR), germline variant calling with HaplotypeCaller, somatic variant calling with Mutect2, and quality metrics collection. The module runs completely standalone with automatic test data download. Individual tasks can be imported and used with existing BAM files and reference data for production workflows.

The module implements GATK best practices for variant calling, including proper duplicate marking, base quality recalibration using known variant sites, **automatic interval splitting for parallelization**, and provides comprehensive validation and reporting for quality assurance.

## Key Features

- **Dual parallelization strategies**: Both scatter-gather and internal parallelization approaches
- **Significant speedup**: Near-linear performance scaling for HaplotypeCaller and Mutect2 on WGS data
- **Smart interval management**: Uses GATK SplitIntervals to create balanced computational chunks
- **Seamless merging**: Automatically combines parallel results into final VCF files
- **Flexible processing modes**: Both individual task execution and combined processing approaches
- **Zero configuration**: Works out-of-the-box with sensible defaults for parallelization

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `mark_duplicates`, `base_recalibrator`, `markdup_recal_metrics`, `haplotype_caller`, `mutect2`, `haplotype_caller_parallel`, `mutect2_parallel`, `split_intervals`, `print_reads`, `merge_vcfs`, `merge_mutect_stats`, `create_sequence_dictionary`, `collect_wgs_metrics`, `fastq_to_bam`, `validate_sam_file`, `analyze_saturation_mutagenesis`, `create_somatic_pon`
- **Test workflow**: `testrun.wdl` (demonstration workflow using test data with parallelization)
- **Container**: `getwilds/gatk:4.6.1.0`

## Tasks

### Core Processing Tasks

### `mark_duplicates`
Marks duplicate reads in aligned BAM files to improve variant calling accuracy.

**Inputs:**
- `bam` (File): Aligned input BAM file
- `bam_index` (File): Index file for the aligned input BAM
- `base_file_name` (String): Base name for the output files
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `markdup_bam` (File): BAM file with duplicate reads marked
- `markdup_bai` (File): Index file for the duplicate-marked BAM
- `duplicate_metrics` (File): Metrics file containing duplicate marking statistics

### `base_recalibrator`
Performs Base Quality Score Recalibration (BQSR) to improve the accuracy of base quality scores.

**Inputs:**
- `bam` (File): Input aligned BAM file to be recalibrated
- `bam_index` (File): Index file for the input BAM
- `dbsnp_vcf` (File): dbSNP VCF file for known variant sites
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `reference_dict` (File): Reference genome sequence dictionary
- `known_indels_sites_vcfs` (Array[File]): Array of VCF files with known indel sites
- `base_file_name` (String): Base name for output files
- `intervals` (File?): Optional interval list file defining target regions
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `recalibrated_bam` (File): BAM file with recalibrated base quality scores
- `recalibrated_bai` (File): Index file for the recalibrated BAM
- `recalibration_report` (File): Base recalibration report table

### `markdup_recal_metrics`
Performs duplicate marking, base recalibration, and WGS metrics collection in a single task to avoid data duplication and improve efficiency.

**Inputs:**
- `bam` (File): Aligned input BAM file
- `bam_index` (File): Index file for the aligned input BAM
- `dbsnp_vcf` (File): dbSNP VCF file for known variant sites
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `reference_dict` (File): Reference genome sequence dictionary
- `known_indels_sites_vcfs` (Array[File]): Array of VCF files with known indel sites
- `base_file_name` (String): Base name for output files
- `intervals` (File?): Optional interval list file defining target regions
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)
- `minimum_mapping_quality` (Int): Minimum mapping quality for reads (default: 20)
- `minimum_base_quality` (Int): Minimum base quality for bases (default: 20)
- `coverage_cap` (Int): Maximum coverage depth to analyze (default: 250)

**Outputs:**
- `recalibrated_bam` (File): BAM file with recalibrated base quality scores
- `recalibrated_bai` (File): Index file for the recalibrated BAM
- `recalibration_report` (File): Base recalibration report table
- `duplicate_metrics` (File): Metrics file containing duplicate marking statistics
- `wgs_metrics` (File): Comprehensive WGS metrics file with coverage and quality statistics

### Core Variant Calling Tasks

### `haplotype_caller`
Calls germline variants using GATK HaplotypeCaller for individual intervals in scatter-gather workflows.

**Inputs:**
- `bam` (File): Input aligned BAM file
- `bam_index` (File): Index file for the input BAM
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `reference_dict` (File): Reference genome sequence dictionary
- `dbsnp_vcf` (File): dbSNP VCF file for variant annotation
- `base_file_name` (String): Base name for output files
- `intervals` (File?): Optional interval list file defining target regions
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `vcf` (File): Compressed VCF file containing germline variant calls
- `vcf_index` (File): Index file for the VCF output

### `haplotype_caller_parallel`
Calls germline variants using GATK HaplotypeCaller with internal parallelization for reduced data duplication and improved efficiency.

**Inputs:**
- `bam` (File): Input aligned BAM file
- `bam_index` (File): Index file for the input BAM
- `intervals` (Array[File]): Array of interval files for parallel processing
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `reference_dict` (File): Reference genome sequence dictionary
- `dbsnp_vcf` (File): dbSNP VCF file for variant annotation
- `base_file_name` (String): Base name for output files
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (should match number of intervals) (default: 2)

**Outputs:**
- `vcf` (File): Compressed VCF file containing germline variant calls
- `vcf_index` (File): Index file for the VCF output

### `mutect2`
Calls somatic variants using GATK Mutect2 in tumor-only mode with filtering for individual intervals in scatter-gather workflows.

**Inputs:**
- `bam` (File): Input aligned BAM file
- `bam_index` (File): Index file for the input BAM
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `reference_dict` (File): Reference genome sequence dictionary
- `gnomad_vcf` (File): gnomAD population allele frequency VCF for germline resource
- `base_file_name` (String): Base name for output files
- `intervals` (File?): Optional interval list file defining target regions
- `max_mnp_distance` (Int): Distance at which to merge MNPs (default: 1, use 0 for PON creation)
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `vcf` (File): Compressed VCF file containing filtered somatic variant calls
- `vcf_index` (File): Index file for the Mutect2 VCF output
- `unfiltered_vcf` (File): Compressed VCF file containing unfiltered somatic variant calls (for PON creation)
- `unfiltered_vcf_index` (File): Index file for the unfiltered Mutect2 VCF output
- `stats_file` (File): Mutect2 statistics file
- `f1r2_counts` (File): F1R2 counts for filtering

### `mutect2_parallel`
Calls somatic variants using GATK Mutect2 with internal parallelization for reduced data duplication and improved efficiency.

**Inputs:**
- `bam` (File): Input aligned BAM file
- `bam_index` (File): Index file for the input BAM
- `intervals` (Array[File]): Array of interval files for parallel processing
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `reference_dict` (File): Reference genome sequence dictionary
- `gnomad_vcf` (File): gnomAD population allele frequency VCF for germline resource
- `base_file_name` (String): Base name for output files
- `max_mnp_distance` (Int): Distance at which to merge MNPs (default: 1, use 0 for PON creation)
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `vcf` (File): Compressed VCF file containing filtered somatic variant calls
- `vcf_index` (File): Index file for the Mutect2 VCF output
- `unfiltered_vcf` (File): Compressed VCF file containing unfiltered somatic variant calls (for PON creation)
- `unfiltered_vcf_index` (File): Index file for the unfiltered Mutect2 VCF output
- `stats_file` (File): Merged Mutect2 statistics file

### `create_somatic_pon`
Creates a somatic panel of normals (PON) from Mutect2 VCF files for filtering germline variants and technical artifacts in tumor-only somatic variant calling.

**Important Requirements:**
- Input VCFs must be **unfiltered** Mutect2 VCFs generated with `max_mnp_distance` set to 0
- VCFs must come from at least two different normal samples

**Inputs:**
- `normal_vcfs` (Array[File]): Array of **unfiltered** Mutect2 VCFs generated with --max-mnp-distance=0
- `normal_vcf_indices` (Array[File]): Array of index files for the normal VCFs
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `reference_dict` (File): Reference genome sequence dictionary
- `intervals` (File): Genomic intervals list, such as from GATK BedToIntervalList
- `database_name` (String): Name for GenomicsDB workspace
- `base_file_name` (String): Base name for output files
- `germline_resource` (File?): Optional gnomAD VCF for additional germline filtering
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `pon_vcf` (File): Gzipped VCF file containing the panel of normals
- `pon_vcf_index` (File): Index file for the panel of normals VCF

### Parallelization and Utility Tasks

### `split_intervals`
Automatically splits genome intervals into optimal chunks for parallel processing using GATK SplitIntervals.

**Inputs:**
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `reference_dict` (File): Reference genome sequence dictionary
- `intervals` (File?): Optional interval list file defining target regions to split
- `scatter_count` (Int): Number of interval files to create (default: 24)
- `filter_to_canonical_chromosomes` (Boolean): Whether to restrict analysis to canonical chromosomes (chr1-22,X,Y,M) (default: true)
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `interval_files` (Array[File]): Array of interval files optimized for parallel processing

### `print_reads`
Extracts reads from specific intervals using GATK PrintReads for scatter-gather processing.

**Inputs:**
- `bam` (File): Input BAM file to extract reads from
- `bam_index` (File): Index file for the input BAM
- `intervals` (Array[File]): Array of interval files defining regions to extract
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `reference_dict` (File): Reference genome sequence dictionary
- `output_basename` (String): Base name for output files
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `interval_bams` (Array[File]): Array of BAM files containing reads from specified intervals
- `interval_bam_indices` (Array[File]): Array of index files for the interval BAMs

### `merge_vcfs`
Merges multiple VCF files from parallel processing into a single consolidated VCF.

**Inputs:**
- `vcfs` (Array[File]): Array of VCF files to merge
- `vcf_indices` (Array[File]): Array of VCF index files
- `base_file_name` (String): Base name for output files
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
- `base_file_name` (String): Base name for output files
- `memory_gb` (Int): Memory allocation in GB (default: 4)
- `cpu_cores` (Int): Number of CPU cores to use (default: 1)

**Outputs:**
- `merged_stats` (File): Merged Mutect2 statistics file

### Utility Tasks

### `fastq_to_bam`
Converts paired FASTQ files to unmapped BAM/SAM format using GATK FastqToSam.

**Inputs:**
- `r1_fastq` (Array[File]): Array of R1 FASTQ files for the library
- `r2_fastq` (Array[File]): Array of R2 FASTQ files for the library
- `base_file_name` (String): Base name for output file
- `sample_name` (String): Sample name to insert into the read group header
- `library_name` (String?): Library name to place into the LB attribute in the read group header (defaults to sample_name if not provided)
- `platform` (String): Sequencing platform (default: illumina)
- `sequencing_center` (String?): Location where the sample was sequenced (defaults to '.' if not provided)
- `read_group_name` (String?): Read group name (if not provided, defaults to sample_name)
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 4)

**Outputs:**
- `unmapped_bam` (File): Unmapped BAM file containing reads from input FASTQ files

**Note:** While `library_name` and `sequencing_center` have defaults for convenience, it's recommended to provide actual values when available for proper read group metadata tracking.

### `validate_sam_file`
Validates BAM/CRAM/SAM files for formatting issues using GATK ValidateSamFile.

**Inputs:**
- `input_file` (File): BAM/CRAM/SAM file to validate
- `base_file_name` (String): Base name for output validation file
- `mode` (String): Validation mode: VERBOSE (detailed), SUMMARY (summary only) (default: SUMMARY)
- `ignore_warnings` (Boolean): Whether to ignore warnings (default: false)
- `reference_fasta` (File?): Reference genome FASTA (required for CRAM files)
- `memory_gb` (Int): Memory allocation in GB (default: 4)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `validation_report` (File): Text file containing validation statistics and any errors/warnings

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
- `base_file_name` (String): Base name for output files
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)
- `minimum_mapping_quality` (Int): Minimum mapping quality for reads (default: 20)
- `minimum_base_quality` (Int): Minimum base quality for bases (default: 20)
- `coverage_cap` (Int): Maximum coverage depth to analyze (default: 250)

**Outputs:**
- `metrics_file` (File): Comprehensive WGS metrics file with coverage and quality statistics

### `analyze_saturation_mutagenesis`
Analyzes saturation mutagenesis data using GATK AnalyzeSaturationMutagenesis to quantify amino acid and codon frequencies.

**Important Requirements:**
- Input BAM will be automatically sorted by queryname (paired reads must be adjacent)
- Reference FASTA must contain **only A, C, G, T bases** (no N's or other IUPAC codes)
- Use `ww-testdata.create_clean_amplicon_reference` to prepare a clean reference if needed

**Inputs:**
- `bam` (File): Input aligned BAM file (will be automatically sorted by queryname)
- `bam_index` (File): Index file for the input BAM
- `reference_fasta` (File): Reference genome FASTA file (must contain only A, C, G, T bases)
- `reference_fasta_index` (File): Index file for the reference FASTA
- `reference_dict` (File): Reference genome sequence dictionary
- `orf_range` (String): Open reading frame range to analyze (e.g., '1-100')
- `base_file_name` (String): Base name for output files
- `memory_gb` (Int): Memory allocation in GB (default: 16)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `aa_counts` (File): Amino acid count table
- `aa_fractions` (File): Amino acid fraction table
- `codon_counts` (File): Codon count table
- `codon_fractions` (File): Codon fraction table
- `cov_length_counts` (File): Coverage length count table
- `read_counts` (File): Read count table
- `ref_coverage` (File): Reference coverage table
- `variant_counts` (File): Variant count table

## Workflow Parallelization Configuration

**Default behavior:**
- Splits genome into configurable intervals (default: 2 for testing, 24 recommended for production)
- Provides both scatter-gather and internal parallelization approaches
- Automatically merges results into final VCF files
- Near-linear speedup for WGS analysis

**Parallelization strategies:**
- **Scatter-gather approach**: Uses `print_reads` to split BAMs by intervals, then runs `haplotype_caller` and `mutect2` on each interval
- **Internal parallelization**: Uses `haplotype_caller_parallel` and `mutect2_parallel` with GNU parallel for more efficient resource usage
- **Combined processing**: `markdup_recal_metrics` performs all preprocessing steps in one task for efficiency

## Testing the Module

The module includes a test workflow that runs with test data and requires no inputs:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint gatk_example
```

### Test Data Workflow

The test workflow automatically:
1. Downloads reference genome data using `ww-testdata`
2. Downloads variant databases (dbSNP, known indels, gnomAD)
3. Downloads test BAM file for demonstration
4. Performs complete GATK variant calling pipeline with parallelization
5. Tests saturation mutagenesis analysis with a small clean reference region
6. Validates all outputs and generates comprehensive reports

No input files or parameters are required - the workflow uses pre-configured test data and parameters optimized for demonstration purposes.

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access for downloading test data (when using automatic mode)
- Sufficient computational resources (memory-intensive for whole genome analysis)
- Multiple CPU cores recommended for optimal parallelization benefits

## Features

- **Dual parallelization strategies**: Both scatter-gather and internal parallelization approaches for maximum flexibility
- **Complete GATK pipeline**: Duplicate marking, base recalibration, germline and somatic variant calling, QC metrics
- **Intelligent interval splitting**: Uses GATK SplitIntervals for optimal load balancing
- **Seamless result merging**: Transparent combination of parallel results
- **Flexible processing modes**: Choose between individual tasks, combined processing, or different parallelization approaches
- **Best practices implementation**: Follows GATK best practices for variant calling workflows
- **Comprehensive validation**: Built-in output validation and quality reporting
- **Resource optimization**: Configurable memory and CPU allocation for each task

## Performance Considerations

### Parallelization Benefits
- **HaplotypeCaller**: Up to 20x speedup on WGS with 24 intervals
- **Mutect2**: Similar dramatic performance improvements
- **Scalability**: Performance scales nearly linearly with available CPU cores
- **Memory efficiency**: Parallel tasks use memory more efficiently than single large jobs
- **Internal parallelization**: The `*_parallel` tasks offer better resource utilization by avoiding data duplication

### Parallelization Strategy Comparison
- **Scatter-gather approach**: Better for distributed computing environments, easier debugging of individual intervals
- **Internal parallelization**: More efficient resource usage, reduced data duplication, better for single-node execution
- **Both approaches**: Available in the same workflow for comparison and flexibility

### Resource Requirements
- **Memory usage**: 16-32GB RAM total for WGS analysis (distributed across parallel tasks)
- **CPU scaling**: **Significant benefit from multiple cores** - recommend 8-24+ CPUs for WGS
- **Storage requirements**: Ensure sufficient disk space for reference data, BAMs, and VCF outputs
- **Network access**: Initial runs require internet connectivity for downloading reference databases
- **Region targeting**: Using intervals files significantly reduces runtime and resource requirements

### Optimization Tips
- **Production workflows**: Import individual tasks and configure scatter_count of 24+ for production data
- **Choose parallelization strategy**: Internal parallelization for single nodes, scatter-gather for distributed computing
- **Ensure adequate CPU allocation**: More cores = better performance for variant calling
- **Monitor memory usage**: Each parallel task needs sufficient memory
- **Choose processing approach**: Combined tasks reduce I/O overhead for large datasets

## Output Description

- **Duplicate-marked BAMs**: BAM files with duplicate reads marked for improved variant calling
- **Recalibrated BAMs**: Quality-improved BAM files with recalibrated base scores
- **Germline VCFs**: HaplotypeCaller variant calls suitable for population genetics and clinical analysis (automatically merged from parallel processing)
- **Somatic VCFs**: Mutect2 tumor-only somatic variant calls with filtering applied (automatically merged from parallel processing)
- **QC metrics**: Comprehensive sequencing quality metrics including coverage statistics and duplicate rates

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Chromosome 1 subset data for efficient testing
- Comprehensive validation of all outputs and statistics
- Integration testing with `ww-testdata` module
- Parallelization testing with multiple interval configurations
- Validation of both scatter-gather and internal parallelization approaches

## Integration Patterns

This module demonstrates several key patterns:
- **Resource management**: Coordinated memory allocation across compute-intensive tasks
- **Best practices workflow**: Implementation of GATK recommended variant calling pipeline
- **Comprehensive validation**: Quality assurance for complex multi-output workflows
- **Dual parallelization approaches**: Both scatter-gather and internal parallelization strategies
- **Flexible processing architectures**: Multiple task orchestration patterns

## Extending the Module

This module can be extended by:
- Adding joint calling capabilities for population-scale analysis
- Integrating additional variant filtering and annotation tools
- Including structural variant calling with other GATK tools
- Adding quality control modules (e.g., FastQC, MultiQC integration)
- Implementing additional parallelization strategies for specific use cases
- Adding performance benchmarking between parallelization approaches

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
