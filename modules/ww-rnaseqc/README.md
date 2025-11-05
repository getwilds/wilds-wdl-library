# ww-rnaseqc Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for running RNA-SeQC quality control analysis on aligned RNA-seq data.

## Overview

This module provides WDL tasks for running RNA-SeQC, a tool for comprehensive quality control metrics on RNA-seq BAM files. RNA-SeQC generates detailed QC metrics including gene coverage, strand specificity, rRNA content, and transcript integrity numbers.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-rnaseqc.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `run_rnaseqc`

Run RNA-SeQC quality control metrics on aligned RNA-seq data.

**Inputs:**
- `bam_file` (File): Aligned reads in BAM format
- `bam_index` (File): Index file for the aligned BAM file
- `ref_gtf` (File): Reference genome GTF annotation file (collapsed)
- `sample_name` (String): Sample name for output files
- `cpu_cores` (Int, default=2): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=4): Memory allocated for the task in GB

**Outputs:**
- `rnaseqc_metrics` (File): Compressed tarball containing RNA-SeQC quality metrics and coverage reports


## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-rnaseqc/ww-rnaseqc.wdl" as rnaseqc_tasks

workflow my_rnaseq_pipeline {
  input {
    File bam_file
    File bam_index
    File gtf_file
    String sample_id
  }

  call rnaseqc_tasks.run_rnaseqc {
    input:
      bam_file = bam_file,
      bam_index = bam_index,
      ref_gtf = gtf_file,
      sample_name = sample_id
  }

  output {
    File qc_metrics = run_rnaseqc.rnaseqc_metrics
  }
}
```

### Advanced Usage Examples

**Custom resource allocation:**
```wdl
call rnaseqc_tasks.run_rnaseqc {
  input:
    bam_file = large_bam,
    bam_index = large_bam_index,
    ref_gtf = reference_gtf,
    sample_name = "high_depth_sample",
    cpu_cores = 4,
    memory_gb = 8
}
```

**Processing multiple samples:**
```wdl
scatter (sample in samples) {
  call rnaseqc_tasks.run_rnaseqc {
    input:
      bam_file = sample.bam,
      bam_index = sample.bai,
      ref_gtf = reference_gtf,
      sample_name = sample.name
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-star**: Use STAR alignment outputs directly as inputs to RNA-SeQC
- **ww-testdata**: Automatic provisioning of test data for demonstrations
- **RNA-seq workflows**: Integrate QC into comprehensive RNA-seq analysis pipelines

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using Sprocket (recommended)
sprocket run testrun.wdl --entrypoint rnaseqc_example

# Using miniWDL
miniwdl run testrun.wdl

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads test reference data using `ww-testdata`
2. Downloads test BAM data using `ww-testdata`
3. Runs RNA-SeQC on the test data
4. Demonstrates the module's tasks in a realistic workflow context

## Docker Container

This module uses the `getwilds/rnaseqc:2.4.2` container image, which includes:
- RNA-SeQC version 2.4.2
- All necessary system dependencies
- Optimized for reproducible RNA-seq QC analysis

## Citation

If you use this module in your research, please cite RNA-SeQC:

> **RNA-SeQC: Comprehensive Quality Control for RNA Sequencing Data**
> DeLuca DS, Levin JZ, Sivachenko A, Fennell T, Nazaire MD, Williams C, Reich M, Winckler W, Getz G.
> Bioinformatics. 2012 Jun 1;28(11):1530-2.
> DOI: [10.1093/bioinformatics/bts196](https://doi.org/10.1093/bioinformatics/bts196)

## RNA-SeQC Output Metrics

RNA-SeQC generates comprehensive quality control metrics including:

### Key Metrics
- **Mapping statistics**: Total reads, mapped reads, unique reads
- **Gene coverage**: 5' to 3' gene body coverage
- **Strand specificity**: Sense vs. antisense mapping rates
- **rRNA content**: Proportion of reads mapping to ribosomal RNA
- **Transcript integrity**: Fragment coverage uniformity
- **Exonic/intronic/intergenic rates**: Read distribution across genomic features

### Output Files
The compressed tarball contains:
- `metrics.tsv`: Tab-delimited summary of all QC metrics
- `coverage.tsv`: Gene body coverage distribution
- Various plots and supplementary files

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 2 cores
- **Memory**: 4 GB
- **Runtime**: ~5-15 minutes per sample depending on BAM size

### Resource Scaling
For larger datasets, consider increasing resources:
- `cpu_cores`: Increase for faster processing (2-4 cores recommended)
- `memory_gb`: Increase for high-depth samples (4-8 GB typically sufficient)

### Input Requirements
- **BAM file**: Must be coordinate-sorted and indexed
- **GTF file**: Must match the reference genome used for alignment
- RNA-SeQC works best with GTF files that have been "collapsed" (merged overlapping transcripts)

## Contributing

To improve this module or report issues:
1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Make your changes following WILDS conventions
3. Test thoroughly with the demonstration workflow
4. Submit a pull request with detailed documentation

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines

## Related Resources

- **[RNA-SeQC GitHub](https://github.com/getzlab/rnaseqc)**: Official RNA-SeQC repository
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
