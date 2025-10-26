# ww-testdata
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

WILDS WDL module for downloading reference and test data.

## Overview

The `ww-testdata` module provides standardized tasks for downloading reference genomes, annotations, and specialized datasets used across WILDS WDL workflows. This module serves as the foundation for consistent, reproducible test data management in the WILDS ecosystem.

## Purpose

Rather than maintaining large static test datasets, `ww-testdata` enables:
- **On-demand data retrieval**: Download only the data you need, when you need it
- **Standardized datasets**: Consistent reference data across all WILDS modules
- **Modular testing**: Each module can specify its exact data requirements
- **Reproducibility**: Versioned, trackable data sources with checksums
- **Efficiency**: Avoid downloading unnecessary test data for targeted module testing

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `download_ref_data`, `download_fastq_data`, `interleave_fastq`, `download_cram_data`, `download_bam_data`, `download_ichor_data`, `download_dbsnp_vcf`, `download_known_indels_vcf`, `download_gnomad_vcf`, `download_annotsv_vcf`, `generate_pasilla_counts`
- **Test workflow**: `testrun.wdl` (demonstration workflow that executes all tasks)

## Usage

### As Imported Tasks

Import specific tasks into your workflows:

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as testdata

workflow my_analysis {
  call testdata.download_ref_data {
    input:
      chromo = "chr22",
      version = "hg38"
  }

  call my_analysis_task {
    input:
      reference_fasta = download_ref_data.fasta,
      reference_gtf = download_ref_data.gtf
  }
}
```

### Test Workflow

### No Input Required

The `testrun.wdl` workflow requires no input parameters and automatically downloads a complete test dataset with hardcoded settings:

- **Chromosome**: chr1 only (for efficient testing)
- **Reference version**: hg38 (latest standard)
- **All test data types**: Reference, FASTQ, CRAM, BAM, ichorCNA, VCF files, and Pasilla counts

### Running the Test Workflow

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Common Integration Patterns

**RNA-seq analysis**:
```wdl
call testdata.download_ref_data {
  input:
    chromo = "chr22",
    version = "hg38"
}
call star_tasks.build_star_index {
  input:
    reference_fasta = download_ref_data.fasta
}
```

**DESeq2 analysis with individual count files**:
```wdl
call testdata.generate_pasilla_counts {
  input:
    n_samples = 6,
    n_genes = 5000
}
call combine_count_matrices {
  input:
    gene_count_files = generate_pasilla_counts.individual_count_files,
    sample_names = generate_pasilla_counts.sample_names,
    sample_conditions = generate_pasilla_counts.sample_conditions
}
call deseq2_analysis {
  input:
    counts_matrix = combine_count_matrices.counts_matrix
}
```

**Copy number analysis**:
```wdl
call testdata.download_ichor_data { }
call ichor_tasks.run_ichor {
  input:
    gc_wig = download_ichor_data.wig_gc,
    map_wig = download_ichor_data.wig_map
}
```

**Variant calling with GATK best practices**:
```wdl
call testdata.download_dbsnp_vcf {
  input:
    region = "chr1:1-10000000",
    filter_name = "chr1"
}
call testdata.download_known_indels_vcf {
  input:
    region = "chr1:1-10000000", 
    filter_name = "chr1"
}
call gatk_tasks.base_recalibrator {
  input:
    known_sites = [download_dbsnp_vcf.dbsnp_vcf, download_known_indels_vcf.known_indels_vcf]
}
```

**Variant annotation**:
```wdl
call testdata.download_annotsv_vcf { }
call annotsv_tasks.annotate_variants {
  input:
    test_vcf = download_annotsv_vcf.test_vcf
}
```

**CRAM processing workflows**:
```wdl
call testdata.download_cram_data {
  input:
    ref_fasta = download_ref_data.fasta
}
call my_cram_analysis {
  input:
    input_cram = download_cram_data.cram,
    input_crai = download_cram_data.crai
}
```

**BAM processing workflows**:
```wdl
call testdata.download_bam_data { }
call my_bam_analysis {
  input:
    input_bam = download_bam_data.bam,
    input_bai = download_bam_data.bai
}
```

## Task Reference

### download_ref_data

**Inputs**:
- `chromo` (String): Chromosome to download
- `version` (String): Genome version
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Note**: In the test workflow, `chromo` is hardcoded to "chr1" and `version` to "hg38".

**Outputs**:
- `fasta` (File): Reference chromosome FASTA file (decompressed)
- `fasta_index` (File): Samtools FASTA index (.fai)
- `gtf` (File): Chromosome-specific gene annotations
- `bed` (File): BED file covering entire chromosome

### download_fastq_data

**Inputs**:
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `r1_fastq` (File): R1 FASTQ file for paired-end sequencing
- `r2_fastq` (File): R2 FASTQ file for paired-end sequencing

### interleave_fastq

**Inputs**:
- `r1_fq` (File): R1 FASTQ file for paired-end sequencing
- `r2_fq` (File): R2 FASTQ file for paired-end sequencing
- `cpu_cores` (Int): CPU allocation (default: 2)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `inter_fastq` (File): Interleaved FASTQ file

### download_cram_data

**Inputs**:
- `ref_fasta` (File): Reference genome FASTA file to use for CRAM conversion
- `cpu_cores` (Int): CPU allocation (default: 2)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `cram` (File): Example CRAM alignment file
- `crai` (File): CRAM index file

### download_bam_data

**Inputs**:
- `filename` (String): Filename to save the BAM file as (default: "NA12878_chr1.bam")
- `cpu_cores` (Int): CPU allocation (default: 2)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `bam` (File): Processed BAM alignment file (chr1 only, primary alignments, 10% subsampled)
- `bai` (File): BAM index file

### download_ichor_data

**Inputs**:
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `wig_gc` (File): GC content in 500kb bins
- `wig_map` (File): Mappability in 500kb bins
- `centromeres` (File): Centromere coordinates
- `panel_of_norm_rds` (File): Panel of normals for normalization

### download_dbsnp_vcf

**Inputs**:
- `region` (String?): Chromosomal region to filter (e.g., "NC_000001.11:1-10000000")
- `filter_name` (String): Filename tag for output (default: "hg38")
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `dbsnp_vcf` (File): Filtered and compressed dbSNP VCF file

### download_known_indels_vcf

**Inputs**:
- `region` (String?): Chromosomal region to filter (e.g., "chr1:1-10000000")
- `filter_name` (String): Filename tag for output (default: "hg38")
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `known_indels_vcf` (File): Filtered and compressed known indels VCF file

### download_gnomad_vcf

**Inputs**:
- `region` (String?): Chromosomal region to filter (e.g., "chr1:1-10000000")
- `filter_name` (String): Filename tag for output (default: "hg38")
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `gnomad_vcf` (File): Filtered and compressed gnomAD VCF file

### download_annotsv_vcf

**Inputs**:
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `test_vcf` (File): Example VCF file for testing

### generate_pasilla_counts

**Inputs**:
- `n_samples` (Int): Number of samples to include (default: 7, max: 7 for pasilla dataset)
- `n_genes` (Int): Approximate number of genes to include (default: 10000)
- `condition_name` (String): Name for the condition column in metadata (default: "condition")
- `output_prefix` (String): Prefix for output files (default: "pasilla")
- `memory_gb` (Int): Memory allocation in GB (default: 4)
- `cpu_cores` (Int): CPU allocation (default: 1)

**Outputs**:
- `individual_count_files` (Array[File]): Individual STAR-format count files for each sample (*.ReadsPerGene.out.tab)
- `sample_names` (Array[String]): Array of sample names corresponding to the count files
- `sample_conditions` (Array[String]): Array of experimental conditions for each sample
- `gene_info` (File): Gene annotation information including gene IDs

## Data Sources

All reference data is downloaded from authoritative public repositories:

- **UCSC Genome Browser**: Reference genomes and annotations
- **GATK Test Data**: Example FASTQ, CRAM, and BAM files
- **ichorCNA Repository**: Copy number analysis references
- **AnnotSV Repository**: Structural variant test data
- **NCBI dbSNP**: Latest dbSNP variant database
- **GATK Resource Bundle**: Known indels and gnomAD population frequencies
- **Bioconductor pasilla package**: Example RNA-seq count data for DESeq2 testing

Data integrity is maintained through the use of stable URLs and version-pinned resources.

## Requirements

### Runtime Dependencies
- **Containers**: `getwilds/samtools:1.11`, `getwilds/awscli:2.27.49`, `getwilds/bcftools:1.19`, `getwilds/deseq2:1.40.2`
- **Tools**: samtools (for FASTA indexing and BAM processing), bcftools (for VCF processing), wget, aws CLI, R with DESeq2 and pasilla packages
- **Network**: Internet access required for data downloads

### Resource Requirements
- **CPU**: 1-2 cores (configurable)
- **Memory**: 2-4 GB (configurable)
- **Storage**: Varies by dataset (chr22: ~50MB, chr1: ~250MB)

## Integration with WILDS Ecosystem

This module is specifically designed to support other WILDS modules:

- **ww-star**: RNA-seq alignment (requires reference FASTA + GTF)
- **ww-bwa**: DNA alignment (requires reference FASTA)
- **ww-deseq2**: Differential expression analysis (uses individual count files from `generate_pasilla_counts`)
- **ww-ichorcna**: Copy number analysis (requires ichorCNA reference files)
- **ww-annotsv**: Structural variant annotation (requires test VCF)
- **Variant calling workflows**: GATK best practices (requires dbSNP, known indels, gnomAD)

By centralizing test data downloads, `ww-testdata` enables:
- Consistent data across all WILDS workflows
- Efficient GitHub Actions testing (download only needed data)
- Simplified module development and testing
- Clear documentation of data dependencies
- Built-in validation to ensure data integrity
- Realistic testing scenarios that mirror production workflows

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Network connectivity validation
- Cross-platform compatibility testing
- Automated output validation

## Support

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
