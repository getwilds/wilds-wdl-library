# ww-rmats-turbo Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for [rMATS-turbo](https://github.com/Xinglab/rmats-turbo) (Replicate Multivariate Analysis of Transcript Splicing), a tool for detecting differential alternative splicing events from RNA-seq data.

## Overview

rMATS-turbo is a fast and accurate computational tool for detecting differential alternative splicing from replicate RNA-seq data. It identifies five types of alternative splicing events:

- **SE**: Skipped exon
- **A5SS**: Alternative 5' splice site
- **A3SS**: Alternative 3' splice site
- **MXE**: Mutually exclusive exons
- **RI**: Retained intron

rMATS-turbo uses a Bayesian framework to calculate the probability that the difference in splicing levels (PSI - Percent Spliced In) between two groups exceeds a user-defined threshold.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-rmats-turbo.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `rmats`

Main task that performs complete differential alternative splicing analysis between two sample groups (combines prep and post steps).

**Inputs:**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `gtf_file` | File | Yes | - | Gene annotation file in GTF format |
| `sample1_bams` | Array[File] | Yes | - | BAM files for sample group 1 (condition 1) |
| `sample2_bams` | Array[File] | Yes | - | BAM files for sample group 2 (condition 2) |
| `read_length` | Int | Yes | - | Length of each read in the RNA-seq data |
| `read_type` | String | No | `paired` | Type of reads: 'paired' or 'single' |
| `library_type` | String | No | `fr-unstranded` | Library type: 'fr-unstranded', 'fr-firststrand', or 'fr-secondstrand' |
| `output_name` | String | No | `rmats_output` | Prefix for output directory name |
| `variable_read_length` | Boolean | No | `false` | Allow reads with different lengths |
| `anchor_length` | Int | No | `1` | Minimum nucleotides mapped to each splice junction end |
| `novel_splice_sites` | Boolean | No | `false` | Enable detection of unannotated splice sites |
| `stat_off` | Boolean | No | `false` | Skip statistical analysis |
| `paired_stats` | Boolean | No | `false` | Use paired statistical model |
| `cstat` | Float | No | `0.0001` | Cutoff splicing difference for null hypothesis test |
| `individual_counts` | Boolean | No | `false` | Output individual count files for each sample |
| `allow_clipping` | Boolean | No | `false` | Allow alignments with soft or hard clipping |
| `min_intron_length` | Int | No | `50` | Minimum intron length for novel splice site detection |
| `max_exon_length` | Int | No | `500` | Maximum exon length for novel splice site detection |
| `cpu_cores` | Int | No | `4` | Number of CPU cores allocated |
| `memory_gb` | Int | No | `16` | Memory allocated in GB |

**Outputs:**
- `output_directory` (File): Tarball containing all rMATS output files including splice event tables and summary statistics

### `rmats_prep`

Preprocesses BAM files and generates .rmats intermediate files. Useful for distributed processing of large datasets.

**Inputs:**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `gtf_file` | File | Yes | - | Gene annotation file in GTF format |
| `sample_bams` | Array[File] | Yes | - | BAM files to preprocess |
| `read_length` | Int | Yes | - | Length of each read |
| `read_type` | String | No | `paired` | Type of reads: 'paired' or 'single' |
| `library_type` | String | No | `fr-unstranded` | Library type |
| `output_name` | String | No | `rmats_prep` | Prefix for output directory name |
| `variable_read_length` | Boolean | No | `false` | Allow reads with different lengths |
| `anchor_length` | Int | No | `1` | Minimum nucleotides mapped to each splice junction end |
| `novel_splice_sites` | Boolean | No | `false` | Enable detection of unannotated splice sites |
| `allow_clipping` | Boolean | No | `false` | Allow alignments with soft or hard clipping |
| `cpu_cores` | Int | No | `4` | Number of CPU cores allocated |
| `memory_gb` | Int | No | `16` | Memory allocated in GB |

**Outputs:**
- `prep_output` (File): Tarball containing .rmats intermediate files

### `rmats_post`

Loads .rmats files and performs alternative splicing event detection and statistical analysis.

**Inputs:**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `gtf_file` | File | Yes | - | Gene annotation file in GTF format |
| `prep_outputs` | Array[File] | Yes | - | Tarballs containing .rmats files from prep step |
| `read_length` | Int | Yes | - | Length of each read |
| `read_type` | String | No | `paired` | Type of reads: 'paired' or 'single' |
| `output_name` | String | No | `rmats_output` | Prefix for output directory name |
| `stat_off` | Boolean | No | `false` | Skip statistical analysis |
| `paired_stats` | Boolean | No | `false` | Use paired statistical model |
| `cstat` | Float | No | `0.0001` | Cutoff splicing difference for null hypothesis test |
| `individual_counts` | Boolean | No | `false` | Output individual count files for each sample |
| `cpu_cores` | Int | No | `4` | Number of CPU cores allocated |
| `memory_gb` | Int | No | `16` | Memory allocated in GB |

**Outputs:**
- `output_directory` (File): Tarball containing all rMATS output files

### `rmats_stat`

Runs statistical analysis on existing rMATS output files.

**Inputs:**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `gtf_file` | File | Yes | - | Gene annotation file in GTF format |
| `existing_output` | File | Yes | - | Tarball containing existing rMATS output files |
| `read_length` | Int | Yes | - | Length of each read |
| `read_type` | String | No | `paired` | Type of reads: 'paired' or 'single' |
| `output_name` | String | No | `rmats_output` | Prefix for output directory name |
| `paired_stats` | Boolean | No | `false` | Use paired statistical model |
| `cstat` | Float | No | `0.0001` | Cutoff splicing difference for null hypothesis test |
| `cpu_cores` | Int | No | `4` | Number of CPU cores allocated |
| `memory_gb` | Int | No | `16` | Memory allocated in GB |

**Outputs:**
- `output_directory` (File): Tarball containing rMATS output files with updated statistical results

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-rmats-turbo/ww-rmats-turbo.wdl" as rmats_tasks

workflow my_splicing_pipeline {
  input {
    File gtf_annotation
    Array[File] control_bams
    Array[File] treatment_bams
    Int read_length
  }

  call rmats_tasks.rmats {
    input:
      gtf_file = gtf_annotation,
      sample1_bams = control_bams,
      sample2_bams = treatment_bams,
      read_length = read_length,
      output_name = "differential_splicing"
  }

  output {
    File splicing_results = rmats.output_directory
  }
}
```

### Advanced Usage Examples

**Strand-specific library with novel splice sites:**
```wdl
call rmats_tasks.rmats {
  input:
    gtf_file = gtf_annotation,
    sample1_bams = control_bams,
    sample2_bams = treatment_bams,
    read_length = 150,
    library_type = "fr-firststrand",
    novel_splice_sites = true,
    output_name = "stranded_analysis"
}
```

**Single-end reads with paired statistical model:**
```wdl
call rmats_tasks.rmats {
  input:
    gtf_file = gtf_annotation,
    sample1_bams = control_bams,
    sample2_bams = treatment_bams,
    read_length = 100,
    read_type = "single",
    paired_stats = true,
    output_name = "single_end_paired_stats"
}
```

**Distributed processing with prep/post workflow:**
```wdl
# Process sample groups separately
scatter (bam in control_bams) {
  call rmats_tasks.rmats_prep as prep_control {
    input:
      gtf_file = gtf_annotation,
      sample_bams = [bam],
      read_length = read_length
  }
}

scatter (bam in treatment_bams) {
  call rmats_tasks.rmats_prep as prep_treatment {
    input:
      gtf_file = gtf_annotation,
      sample_bams = [bam],
      read_length = read_length
  }
}

# Combine and analyze
call rmats_tasks.rmats_post {
  input:
    gtf_file = gtf_annotation,
    prep_outputs = flatten([prep_control.prep_output, prep_treatment.prep_output]),
    read_length = read_length
}
```

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using Sprocket
sprocket run testrun.wdl --entrypoint rmats_turbo_example

# Using miniWDL
miniwdl run testrun.wdl

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads reference GTF annotation using `ww-testdata`
2. Downloads test BAM files from the GATK test data bucket
3. Runs rMATS-turbo to demonstrate module functionality
4. Tests the prep and stat tasks as well

**Note:** The test BAM files are WGS data, not RNA-seq, so the splicing analysis results will not be biologically meaningful but serve to validate that the module executes correctly.

## Docker Container

This module uses the `getwilds/rmats-turbo:latest` container image from the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library), which includes:
- rMATS-turbo
- Python 3.9+
- All required dependencies

## Input Requirements

### BAM Files

- BAM files should be aligned using a splice-aware aligner (e.g., STAR, HISAT2)
- Files should be sorted and ideally indexed (.bai)
- The same reference genome used for alignment should match the GTF annotation

### GTF Annotation

- Gene annotation in GTF format (e.g., from Ensembl or GENCODE)
- Should match the reference genome used for BAM alignment
- Exon features are required for splice site detection

## Output Files

The output tarball contains:

| File | Description |
|------|-------------|
| `SE.MATS.JC.txt` | Skipped exon events (junction counts only) |
| `SE.MATS.JCEC.txt` | Skipped exon events (junction + exon counts) |
| `A5SS.MATS.JC.txt` | Alternative 5' splice site events (junction counts only) |
| `A5SS.MATS.JCEC.txt` | Alternative 5' splice site events (junction + exon counts) |
| `A3SS.MATS.JC.txt` | Alternative 3' splice site events (junction counts only) |
| `A3SS.MATS.JCEC.txt` | Alternative 3' splice site events (junction + exon counts) |
| `MXE.MATS.JC.txt` | Mutually exclusive exon events (junction counts only) |
| `MXE.MATS.JCEC.txt` | Mutually exclusive exon events (junction + exon counts) |
| `RI.MATS.JC.txt` | Retained intron events (junction counts only) |
| `RI.MATS.JCEC.txt` | Retained intron events (junction + exon counts) |
| `fromGTF.*.txt` | Events identified from GTF annotation |
| `summary.txt` | Summary statistics of detected events |

## Citation

If you use rMATS-turbo in your research, please cite:

> Shen S, Park JW, Lu ZX, et al. rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. Proc Natl Acad Sci U S A. 2014;111(51):E5593-E5601. doi:10.1073/pnas.1419161111

## About rMATS-turbo

- **Purpose**: Detect differential alternative splicing between conditions
- **Input**: Aligned RNA-seq BAM files and GTF annotation
- **Output**: Tables of alternative splicing events with statistical significance
- **Documentation**: https://github.com/Xinglab/rmats-turbo
- **Source Code**: https://github.com/Xinglab/rmats-turbo

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines

## Related Resources

- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
- **[rMATS Documentation](https://github.com/Xinglab/rmats-turbo)**: Official rMATS-turbo documentation
- **[ww-jcast](../ww-jcast/)**: JCAST module for translating rMATS output to protein sequences
