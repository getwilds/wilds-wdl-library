# ww-rseqc Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for running comprehensive RNA-seq quality control analysis using RSeQC.

## Overview

This module provides WDL tasks for running RSeQC, a comprehensive toolkit for RNA-seq quality control. RSeQC generates detailed QC metrics including read distribution, gene body coverage, strand specificity, and splice junction analysis.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-rseqc.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `run_rseqc`

Run comprehensive RSeQC quality control metrics on aligned RNA-seq data.

**Inputs:**
- `bam_file` (File): Aligned reads in BAM format
- `bam_index` (File): Index file for the aligned BAM file
- `ref_bed` (File): Reference genome annotation in BED format (12-column)
- `sample_name` (String): Sample name for output files
- `cpu_cores` (Int, default=2): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=4): Memory allocated for the task in GB

**Outputs:**
- `read_distribution` (File): Distribution of reads across genomic features
- `gene_body_coverage` (File): Gene body coverage metrics (text file)
- `gene_body_coverage_plot` (File): Gene body coverage plot (PDF)
- `infer_experiment` (File): Strand specificity inference results
- `bam_stat` (File): Basic BAM alignment statistics
- `junction_xls` (File): Splice junction annotation in XLS format
- `junction_log` (File): Junction annotation log
- `rseqc_summary` (File): Summary report of all RSeQC metrics

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-rseqc/ww-rseqc.wdl" as rseqc_tasks

workflow my_rnaseq_pipeline {
  input {
    File bam_file
    File bam_index
    File bed_file
    String sample_id
  }

  call rseqc_tasks.run_rseqc {
    input:
      bam_file = bam_file,
      bam_index = bam_index,
      ref_bed = bed_file,
      sample_name = sample_id
  }

  output {
    File read_dist = run_rseqc.read_distribution
    File coverage = run_rseqc.gene_body_coverage
    File strand_info = run_rseqc.infer_experiment
    File qc_summary = run_rseqc.rseqc_summary
  }
}
```

### Advanced Usage Examples

**Custom resource allocation:**
```wdl
call rseqc_tasks.run_rseqc {
  input:
    bam_file = large_bam,
    bam_index = large_bam_index,
    ref_bed = annotation_bed,
    sample_name = "high_depth_sample",
    cpu_cores = 4,
    memory_gb = 8
}
```

**Processing multiple samples:**
```wdl
scatter (sample in samples) {
  call rseqc_tasks.run_rseqc {
    input:
      bam_file = sample.bam,
      bam_index = sample.bai,
      ref_bed = reference_bed,
      sample_name = sample.name
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-star**: Use STAR alignment outputs directly as inputs to RSeQC
- **ww-testdata**: Automatic provisioning of test data for demonstrations
- **RNA-seq workflows**: Integrate QC into comprehensive RNA-seq analysis pipelines

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using Sprocket (recommended)
sprocket run testrun.wdl --entrypoint rseqc_example

# Using miniWDL
miniwdl run testrun.wdl

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads test reference data using `ww-testdata`
2. Downloads test BAM data using `ww-testdata`
3. Converts GTF to BED format for RSeQC
4. Runs RSeQC on the test data
5. Demonstrates the module's tasks in a realistic workflow context

## Docker Container

This module uses the **`getwilds/rseqc:5.0.4`** container image from the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library), which includes:
- RSeQC version 5.0.4
- Python 3.9 environment
- All necessary dependencies
- Optimized for reproducible RNA-seq QC analysis

## Citation

If you use this module in your research, please cite RSeQC:

> **RSeQC: quality control of RNA-seq experiments**
> Wang L, Wang S, Li W.
> Bioinformatics. 2012 Aug 15;28(16):2184-5.
> DOI: [10.1093/bioinformatics/bts356](https://doi.org/10.1093/bioinformatics/bts356)

## RSeQC Analyses Performed

### 1. Read Distribution
Analyzes how reads are distributed across different genomic features:
- CDS (coding sequence)
- 5' UTR and 3' UTR
- Introns
- Intergenic regions

### 2. Gene Body Coverage
Assesses the uniformity of coverage along gene bodies:
- Detects 5' or 3' coverage bias
- Important for assessing RNA degradation
- Generates coverage curve plots

### 3. Infer Experiment
Determines the strand-specificity of the RNA-seq library:
- Unstranded
- Sense-stranded
- Antisense-stranded

### 4. BAM Statistics
Provides basic alignment statistics:
- Total reads
- Mapped reads
- Properly paired reads
- Duplicate reads

### 5. Junction Annotation
Classifies and analyzes splice junctions:
- Known splice junctions
- Novel splice junctions
- Partial novel junctions
- Junction plots and statistics

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
- **BED file**: Must be in 12-column BED format
  - Can be generated from GTF using tools like `gtf2bed` or UCSC utilities
  - Must match the reference genome used for alignment

## Converting GTF to BED

RSeQC requires BED12 format. You can convert GTF to BED using various tools:

```bash
# Using UCSC gtfToGenePred and genePredToBed
gtfToGenePred input.gtf output.genePred
genePredToBed output.genePred output.bed

# Using bedops
gtf2bed < input.gtf > output.bed
```

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

- **[RSeQC Documentation](http://rseqc.sourceforge.net/)**: Official RSeQC documentation and tutorials
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
