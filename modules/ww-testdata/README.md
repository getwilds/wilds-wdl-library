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

- **Tasks**: `download_ref_data`, `download_fastq_data`, `download_test_transcriptome`, `interleave_fastq`, `download_cram_data`, `download_bam_data`, `download_ichor_data`, `download_dbsnp_vcf`, `download_known_indels_vcf`, `download_gnomad_vcf`, `download_annotsv_vcf`, `generate_pasilla_counts`, `create_clean_amplicon_reference`, `create_gdc_manifest`
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
      version = "hg38",
      region = "1-20000000"  # Optional: limit to first 20Mb for faster downloads
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
- **All test data types**: Reference genome, transcriptome, FASTQ, interleaved FASTQ, CRAM, BAM, ichorCNA files, VCF files (dbSNP, known indels, gnomAD, AnnotSV), and Pasilla counts

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

**RNA-seq alignment analysis**:
```wdl
call testdata.download_ref_data {
  input:
    chromo = "chr22",
    version = "hg38",
    region = "1-30000000"  # Optional: subset for faster testing
}
call star_tasks.build_star_index {
  input:
    reference_fasta = download_ref_data.fasta,
    reference_gtf = download_ref_data.gtf
}
```

**RNA-seq quantification analysis**:
```wdl
call testdata.download_test_transcriptome { }
call testdata.download_fastq_data { }
call salmon_tasks.build_index {
  input:
    transcriptome_fasta = download_test_transcriptome.transcriptome_fasta
}
call salmon_tasks.quantify {
  input:
    salmon_index_dir = build_index.salmon_index,
    fastq_r1 = download_fastq_data.r1_fastq,
    fastq_r2 = download_fastq_data.r2_fastq
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
- `chromo` (String): Chromosome to download (default: "chr1")
- `version` (String): Genome version (default: "hg38")
- `region` (String, optional): Region coordinates to extract from chromosome in format '1-30000000'. If not specified, uses entire chromosome
- `output_name` (String, optional): Name for output files (default: uses chromo name)
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Note**: In the test workflow, `chromo` is hardcoded to "chr1", `version` to "hg38", and `region` to "1-10000000" for faster testing.

**Outputs**:
- `fasta` (File): Reference chromosome FASTA file (decompressed, filtered to region if specified)
- `fasta_index` (File): Samtools FASTA index (.fai)
- `dict` (File): Samtools FASTA dictionary file (.dict) for GATK compatibility
- `gtf` (File): Chromosome-specific gene annotations (filtered to region if specified)
- `bed` (File): BED file covering the entire chromosome or specified region

### download_fastq_data

**Inputs**:
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `r1_fastq` (File): R1 FASTQ file for paired-end sequencing
- `r2_fastq` (File): R2 FASTQ file for paired-end sequencing

### download_test_transcriptome

Downloads protein-coding transcriptome from GENCODE for RNA-seq quantification testing (e.g., Salmon, Kallisto).

**Important Note**: This task uses GENCODE (Ensembl) annotations, while other ww-testdata tasks (`download_ref_data`) use NCBI RefSeq annotations. These annotation sources have different gene/transcript IDs and may differ in transcript models. For testing purposes, this provides functional validation of quantification tools. For production pipelines, ensure you maintain annotation consistency throughout your workflow (i.e., use the same annotation source for alignment, quantification, and downstream analysis).

**Inputs**:
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 2)

**Outputs**:
- `transcriptome_fasta` (File): Protein-coding transcriptome FASTA file (~150MB uncompressed, ~20,000 transcripts from GENCODE release 47)

### create_gdc_manifest

Creates a test GDC manifest file containing small open-access files for testing the ww-gdc module. This task generates a properly formatted tab-separated manifest file with file UUIDs, filenames, MD5 checksums, file sizes, and release status.

**Inputs**: None

**Outputs**:
- `manifest` (File): GDC manifest file containing 3 small open-access TCGA files (total ~210KB)

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

### create_clean_amplicon_reference

Extracts and cleans a reference sequence region for saturation mutagenesis analysis. This task is designed to prepare reference sequences for tools that require only standard nucleotides (A, C, G, T), such as GATK AnalyzeSaturationMutagenesis.

**Use Case**: When performing saturation mutagenesis or deep mutational scanning experiments on specific genomic regions (amplicons), tools like GATK's AnalyzeSaturationMutagenesis require reference sequences without ambiguous bases (N's or other IUPAC codes). This task extracts your target region and ensures it contains only A, C, G, T bases.

**Inputs**:
- `input_fasta` (File): Input reference FASTA file
- `region` (String?): Region to extract in format 'chr:start-end' (e.g., 'chr1:1000-2000'). If not specified, uses entire sequence.
- `output_name` (String): Name for the output reference (default: "amplicon")
- `replace_n_with` (String): Base to replace N's with (default: "A"). Use empty string to fail if N's are found.
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 2)

**Outputs**:
- `clean_fasta` (File): Cleaned reference FASTA file with no ambiguous bases
- `clean_fasta_index` (File): Samtools index for the cleaned reference
- `clean_dict` (File): Samtools dictionary for the cleaned reference

**Example Usage**:
```wdl
# For saturation mutagenesis on a specific amplicon
call testdata.create_clean_amplicon_reference {
  input:
    input_fasta = "hg38_chr1.fa",
    region = "chr1:12345-67890",  # Your amplicon coordinates
    output_name = "my_amplicon",
    replace_n_with = "A"  # Replace any N's with A
}

call gatk.analyze_saturation_mutagenesis {
  input:
    reference_fasta = create_clean_amplicon_reference.clean_fasta,
    reference_fasta_index = create_clean_amplicon_reference.clean_fasta_index,
    reference_dict = create_clean_amplicon_reference.clean_dict
}
```

## Data Sources

All reference data is downloaded from authoritative public repositories:

- **UCSC Genome Browser**: Reference genomes and annotations
- **GENCODE**: Human transcriptome annotations and sequences
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
- **Tools**: samtools (for FASTA indexing and BAM processing), bcftools (for VCF processing), curl, aws CLI, R with DESeq2 and pasilla packages
- **Network**: Internet access required for data downloads

### Resource Requirements
- **CPU**: 1-2 cores (configurable)
- **Memory**: 2-4 GB (configurable)
- **Storage**: Varies by dataset (chr22: ~50MB, chr1: ~250MB)

## Integration with WILDS Ecosystem

This module is specifically designed to support other WILDS modules:

- **ww-star**: RNA-seq alignment (requires reference FASTA + GTF)
- **ww-salmon**: RNA-seq quantification (requires transcriptome FASTA)
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
