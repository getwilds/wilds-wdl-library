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

- **Tasks**: `download_ref_data`, `download_fastq_data`, `download_test_transcriptome`, `interleave_fastq`, `download_cram_data`, `download_bam_data`, `download_ichor_data`, `download_dbsnp_vcf`, `download_known_indels_vcf`, `download_gnomad_vcf`, `download_annotsv_vcf`, `generate_pasilla_counts`, `create_clean_amplicon_reference`, `create_gdc_manifest`, `download_shapemapper_data`, `download_test_cellranger_ref`, `download_diamond_data`, `download_glimpse2_genetic_map`, `download_glimpse2_reference_panel`, `download_glimpse2_test_gl_vcf`
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
- **All test data types**: Reference genome, transcriptome, FASTQ, interleaved FASTQ, CRAM, BAM, ichorCNA files, VCF files (dbSNP, known indels, gnomAD, AnnotSV), Pasilla counts, ShapeMapper data, Cell Ranger reference, DIAMOND data, and GLIMPSE2 imputation data

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

**RNA structure analysis with ShapeMapper**:
```wdl
call testdata.download_shapemapper_data { }
call shapemapper_tasks.run_shapemapper {
  input:
    sample_name = "TPP_riboswitch",
    target_fa = download_shapemapper_data.target_fa,
    modified_r1 = download_shapemapper_data.modified_r1,
    modified_r2 = download_shapemapper_data.modified_r2,
    untreated_r1 = download_shapemapper_data.untreated_r1,
    untreated_r2 = download_shapemapper_data.untreated_r2,
    is_amplicon = true,
    min_depth = 1000
}
```

**Single-cell RNA-seq with Cell Ranger**:
```wdl
call testdata.download_test_cellranger_ref { }
call cellranger_tasks.run_count {
  input:
    r1_fastqs = my_r1_fastqs,
    r2_fastqs = my_r2_fastqs,
    ref_gex = download_test_cellranger_ref.ref_tar,
    sample_id = "my_sample"
}
```

**Protein sequence alignment with DIAMOND**:
```wdl
call testdata.download_diamond_data { }
call diamond_tasks.make_database {
  input:
    fasta = download_diamond_data.reference
}
call diamond_tasks.diamond_blastp {
  input:
    db = make_database.db,
    query = download_diamond_data.query
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

Downloads small example FASTQ files for testing. Renames to Illumina naming convention with optional gzip compression.

**Inputs**:
- `prefix` (String): Sample prefix for output filenames (default: "testdata")
- `gzip_output` (Boolean): Compress output files with gzip (default: false)
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `r1_fastq` (File): R1 FASTQ file named `<prefix>_S1_L001_R1_001.fastq[.gz]`
- `r2_fastq` (File): R2 FASTQ file named `<prefix>_S1_L001_R2_001.fastq[.gz]`

**Example Usage**:
```wdl
call testdata.download_fastq_data { input:
  prefix = "my_sample",
  gzip_output = true
}
# Outputs: my_sample_S1_L001_R1_001.fastq.gz, my_sample_S1_L001_R2_001.fastq.gz
```

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

### download_shapemapper_data

Downloads the official ShapeMapper example data (TPP riboswitch) from the Weeks-UNC/shapemapper2 repository. This dataset contains RNA structure probing data suitable for testing ShapeMapper analysis workflows.

**Use Case**: When testing RNA structure analysis workflows with ShapeMapper, you need appropriate test data with actual chemical probing experiments. This task downloads the canonical TPP riboswitch example data that includes both modified (chemically treated) and untreated control samples.

**Inputs**:
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `target_fa` (File): Target RNA FASTA file (TPP riboswitch sequence, ~200 nucleotides)
- `modified_r1` (File): R1 FASTQ file from modified/treated sample (TPPplus, concatenated from split files)
- `modified_r2` (File): R2 FASTQ file from modified/treated sample (TPPplus, concatenated from split files)
- `untreated_r1` (File): R1 FASTQ file from untreated control sample (TPPminus, concatenated from split files)
- `untreated_r2` (File): R2 FASTQ file from untreated control sample (TPPminus, concatenated from split files)

**Data Source**: https://github.com/Weeks-UNC/shapemapper2/tree/master/example_data

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
- `test_vcf` (File): Example VCF file for testing structural variant annotation

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

### download_test_cellranger_ref

Downloads a minimal Cell Ranger reference transcriptome for testing single-cell RNA-seq workflows. This reference contains only chromosomes 21 and 22, making it small enough for CI testing while still being functional.

**Use Case**: When testing Cell Ranger workflows, you need a properly formatted reference transcriptome. Full references are very large (>10GB), but this minimal reference (~728MB) is sufficient for validating workflow execution without excessive download times or storage requirements.

**Inputs**:
- `cpu_cores` (Int): CPU allocation (default: 2)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `ref_tar` (File): Cell Ranger reference transcriptome tarball containing chromosomes 21 and 22

**Data Source**: Swiss Institute of Bioinformatics single-cell training materials (https://sib-swiss.github.io/single-cell-training-archived/)

**Example Usage**:
```wdl
call testdata.download_test_cellranger_ref { }
call cellranger_tasks.run_count {
  input:
    r1_fastqs = my_r1_fastqs,
    r2_fastqs = my_r2_fastqs,
    ref_gex = download_test_cellranger_ref.ref_tar,
    sample_id = "my_sample"
}
```

### download_diamond_data

Downloads E. coli Swiss-Prot reference proteome and creates a small subset for testing DIAMOND protein alignment workflows.

**Use Case**: When testing DIAMOND protein sequence alignment workflows, you need both a reference database and query sequences. This task downloads the E. coli reference proteome from UniProt and creates a small subset of 10 sequences that can be used as a test query.

**Inputs**:
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 2)

**Outputs**:
- `reference` (File): Full E. coli reference proteome FASTA file (ecoli_proteins.fasta, ~4,400 protein sequences)
- `query` (File): Subset of first 10 sequences (ecoli_subset.fasta) for use as test queries

**Example Usage**:
```wdl
# For testing DIAMOND protein alignment
call testdata.download_diamond_data { }
call diamond_tasks.make_database {
  input:
    fasta = download_diamond_data.reference
}
call diamond_tasks.diamond_blastp {
  input:
    db = make_database.db,
    query = download_diamond_data.query,
    align_id = "50",
    query_cover = "50"
}
```

### download_glimpse2_genetic_map

Downloads genetic map files for GLIMPSE2 imputation from the official GLIMPSE repository.

**Use Case**: GLIMPSE2 requires genetic map files for accurate imputation. This task downloads chromosome-specific genetic maps that define recombination rates across the genome.

**Inputs**:
- `chromosome` (String): Chromosome to download genetic map for (default: "chr1")
- `genome_build` (String): Genome build version, "b37" or "b38" (default: "b38")
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 2)

**Outputs**:
- `genetic_map` (File): Compressed genetic map file for the specified chromosome

**Data Source**: https://github.com/odelaneau/GLIMPSE/tree/master/maps

### download_glimpse2_reference_panel

Downloads and prepares a 1000 Genomes reference panel subset for GLIMPSE2 testing. This task downloads phased data for a specified chromosome and filters to a region for efficient CI/CD testing.

**Use Case**: GLIMPSE2 imputation requires a phased reference panel. This task downloads the 1000 Genomes high-coverage phased panel, filters to biallelic SNPs, and creates a sites-only VCF for genotype likelihood calculation.

**Inputs**:
- `chromosome` (String): Chromosome to download (default: "chr1")
- `region` (String): Genomic region to extract (default: "chr1:1-10000000"). Must match the chromosome parameter.
- `exclude_samples` (String): Comma-separated list of samples to exclude, useful for leave-one-out validation (default: "NA12878")
- `cpu_cores` (Int): CPU allocation (default: 2)
- `memory_gb` (Int): Memory allocation (default: 8)

**Outputs**:
- `reference_vcf` (File): Reference panel BCF file for imputation
- `reference_vcf_index` (File): Index file for the reference panel
- `sites_vcf` (File): Sites-only VCF for genotype likelihood calculation
- `sites_vcf_index` (File): Index file for sites VCF

**Data Source**: 1000 Genomes high-coverage phased data (http://ftp.1000genomes.ebi.ac.uk/)

### download_glimpse2_test_gl_vcf

Downloads low-coverage sequencing data from 1000 Genomes and extracts a VCF with genotype likelihoods for GLIMPSE2 imputation testing. Uses NA12878 data from the 1000 Genomes Phase 3 low-coverage dataset.

**Use Case**: GLIMPSE2 can impute from VCF files containing genotype likelihoods (GL fields). This task downloads pre-computed genotype likelihoods from 1000 Genomes rather than generating them from BAM data, providing a simpler and more reliable test data source.

**Inputs**:
- `chromosome` (String): Chromosome to download (default: "chr1"). Note: Phase 3 data uses numeric chromosome names (1-22).
- `region` (String): Genomic region to extract (default: "chr1:1-10000000"). Must match the chromosome parameter.
- `sample_name` (String): Sample to extract from 1000 Genomes (default: "NA12878")
- `cpu_cores` (Int): CPU allocation (default: 2)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `gl_vcf` (File): VCF file with genotype likelihoods (GL field) for imputation
- `gl_vcf_index` (File): Index file for the GL VCF

**Data Source**: 1000 Genomes Phase 3 low-coverage data (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)

**Example Usage**:
```wdl
# For testing GLIMPSE2 imputation on chr22
call testdata.download_glimpse2_test_gl_vcf {
  input:
    chromosome = "chr22",
    region = "chr22:20000000-21000000"
}

call glimpse2_tasks.glimpse2_phase {
  input:
    input_vcf = download_glimpse2_test_gl_vcf.gl_vcf,
    reference_chunk = my_reference_chunk
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
- **EBI 1000 Genomes FTP**: Mills and 1000G gold standard known indels
- **GATK Resource Bundle**: gnomAD population frequencies
- **GLIMPSE Repository**: Genetic maps for imputation (https://github.com/odelaneau/GLIMPSE)
- **1000 Genomes High-Coverage**: Phased reference panels for GLIMPSE2 imputation
- **1000 Genomes Phase 3**: Low-coverage sequencing data with genotype likelihoods for imputation testing
- **Bioconductor pasilla package**: Example RNA-seq count data for DESeq2 testing
- **ShapeMapper Repository**: TPP riboswitch RNA structure probing example data
- **Swiss Institute of Bioinformatics**: Minimal Cell Ranger reference (chr21/22) for single-cell testing
- **UniProt**: E. coli K-12 reference proteome for DIAMOND protein alignment testing

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
- **ww-shapemapper**: RNA structure analysis (uses TPP riboswitch example data from `download_shapemapper_data`)
- **ww-cellranger**: Single-cell RNA-seq analysis (uses minimal reference from `download_test_cellranger_ref`)
- **ww-diamond**: Protein sequence alignment (uses E. coli proteome from `download_diamond_data`)
- **ww-annovar**: Variant annotation (uses gnomAD VCF from `download_gnomad_vcf`)
- **ww-glimpse2**: Genotype imputation (uses genetic maps, reference panels, and GL VCFs from GLIMPSE2 tasks)
- **ww-consensus**: Consensus variant calling (uses gnomAD VCF from `download_gnomad_vcf`)
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
