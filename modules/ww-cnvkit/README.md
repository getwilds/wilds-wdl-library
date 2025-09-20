# ww-cnvkit Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for copy number variation (CNV) detection using CNVkit. This module implements CNVkit workflows for both germline and somatic CNV calling from targeted DNA sequencing data.

## Overview

CNVkit is a Python library and command-line software toolkit to infer and visualize copy number from high-throughput DNA sequencing data. This module provides a comprehensive WDL wrapper for CNVkit functionality, supporting both tumor-only and paired tumor/normal analyses.

**Key Features:**
- Reference creation from normal samples or pooled controls
- Tumor-only and paired tumor/normal CNV analysis
- Automated target and antitarget region handling
- Comprehensive visualization output
- Integration with WILDS ecosystem

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-cnvkit.wdl` - Contains all tasks and demonstration workflow
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `create_reference`

Creates a CNVkit reference from normal samples or pooled reference data.

**Inputs:**
- `bam_files` (Array[File]): Array of BAM files for reference creation
- `bam_indices` (Array[File]): Array of BAM index files corresponding to BAM files
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Reference genome FASTA index file
- `target_bed` (File, optional): Target regions BED file
- `antitarget_bed` (File, optional): Antitarget regions BED file
- `cpu_cores` (Int, default=4): Number of CPU cores to use
- `memory_gb` (Int, default=16): Memory allocation in GB

**Outputs:**
- `reference_cnn` (File): CNVkit reference file (.cnn)

### `run_cnvkit`

Performs CNVkit copy number analysis on tumor samples.

**Inputs:**
- `sample_name` (String): Sample identifier
- `tumor_bam` (File): Tumor BAM file
- `tumor_bai` (File): Tumor BAM index file
- `normal_bam` (File, optional): Normal/control BAM file
- `normal_bai` (File, optional): Normal/control BAM index file
- `reference_cnn` (File): CNVkit reference file
- `target_bed` (File, optional): Target regions BED file
- `paired_analysis` (Boolean, default=false): Whether to perform paired analysis
- `cpu_cores` (Int, default=4): Number of CPU cores to use
- `memory_gb` (Int, default=16): Memory allocation in GB

**Outputs:**
- `cnv_segments` (File): CNV segments file (.cnr)
- `cnv_calls` (File): CNV calls file (.cns)
- `cnv_plot` (File): CNV visualization plot (PDF)

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-cnvkit/ww-cnvkit.wdl" as cnvkit

workflow my_cnv_analysis {
  input {
    Array[File] normal_bams
    Array[File] normal_bais
    File tumor_bam
    File tumor_bai
    File reference_fasta
    File reference_fasta_index
  }

  # Create reference from normal samples
  call cnvkit.create_reference {
    input:
      bam_files = normal_bams,
      bam_indices = normal_bais,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index
  }

  # Analyze tumor sample
  call cnvkit.run_cnvkit {
    input:
      sample_name = "tumor_sample",
      tumor_bam = tumor_bam,
      tumor_bai = tumor_bai,
      reference_cnn = create_reference.reference_cnn
  }

  output {
    File cnv_calls = run_cnvkit.cnv_calls
    File cnv_segments = run_cnvkit.cnv_segments
    File cnv_plot = run_cnvkit.cnv_plot
    File reference_file = create_reference.reference_cnn
  }
}
```

### Advanced Usage Examples

**Tumor-only analysis:**
```wdl
call cnvkit.run_cnvkit {
  input:
    sample_name = "tumor_sample_01",
    tumor_bam = tumor_bam_file,
    tumor_bai = tumor_bai_file,
    reference_cnn = reference_file,
    paired_analysis = false
}
```

**Paired tumor/normal analysis:**
```wdl
call cnvkit.run_cnvkit {
  input:
    sample_name = "paired_sample_01",
    tumor_bam = tumor_bam_file,
    tumor_bai = tumor_bai_file,
    normal_bam = normal_bam_file,
    normal_bai = normal_bai_file,
    reference_cnn = reference_file,
    paired_analysis = true
}
```

**Custom resource allocation:**
```wdl
call cnvkit.run_cnvkit {
  input:
    sample_name = "large_sample",
    tumor_bam = large_tumor_bam,
    tumor_bai = large_tumor_bai,
    reference_cnn = reference_file,
    cpu_cores = 8,
    memory_gb = 32
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-testdata**: Automatic provisioning of test BAM data for demonstrations
- **ww-gatk**: Can be used upstream for BAM preprocessing
- **ww-samtools**: Compatible with samtools-processed BAM files

## Testing the Module

The module includes a demonstration workflow that can be tested independently:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-cnvkit.wdl

# Using miniWDL
miniwdl run ww-cnvkit.wdl

# Using Sprocket
sprocket run ww-cnvkit.wdl
```

### Automatic Demo Mode

The `cnvkit_example` workflow automatically:
1. Downloads reference genome data (FASTA and index) using `ww-testdata`
2. Downloads multiple test BAM files with unique sample names for reference creation
3. Creates a CNVkit reference from the normal samples using WGS mode
4. Processes a single tumor sample against the reference
5. Generates CNV calls (.cns), segments (.cnr), and visualization plots (PDF)

## Docker Container

This module uses the `getwilds/cnvkit:0.9.10` container image, which includes:
- CNVkit version 0.9.10
- Python 3.10 with scientific computing stack (NumPy, SciPy, Pandas, Matplotlib)
- R with essential Bioconductor packages (DNAcopy, PSCBS, GenomicRanges)
- Additional R packages for enhanced plotting (ggplot2, gridExtra, RColorBrewer)
- All necessary system dependencies

## CNVkit Information

- **Tool**: CNVkit copy number analysis toolkit
- **Version**: 0.9.10
- **Purpose**: Copy number variation detection from targeted sequencing
- **Input Types**: BAM files from targeted DNA sequencing
- **Output Types**: CNV segments (.cnr), CNV calls (.cns), visualization plots

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 4 cores
- **Memory**: 16 GB
- **Runtime**: 10-30 minutes per sample depending on data size

### Resource Scaling
CNVkit resource requirements depend on:
- `cpu_cores`: More cores improve parallel processing performance
- `memory_gb`: Increase for large BAM files or whole-genome data
- **BAM file size**: Larger files require more memory and time
- **Target region size**: More regions increase computational complexity

### Performance Guidelines
- **Exome data**: Default resources (4 cores, 16 GB) are typically sufficient
- **Large gene panels**: May benefit from increased memory (24-32 GB)
- **Whole genome**: Requires substantial resources (8+ cores, 32+ GB memory)

## CNVkit Analysis Types

### Reference Creation
- **Pooled normal**: Uses multiple normal samples to create robust reference
- **Process-matched**: Uses samples processed with identical protocols
- **Auto-target**: Automatically identifies target regions from BAM files

### CNV Detection Methods
- **Threshold-based**: Default method using statistical thresholds (currently implemented)

*Note: CNVkit supports additional segmentation methods like CBS (Circular Binary Segmentation) and HMM (Hidden Markov Model), but these are not yet exposed as configurable options in this WDL module.*

## Output Files

### CNV Segments (.cnr)
- Raw copy number ratios for each genomic bin
- Used for visualization and downstream analysis

### CNV Calls (.cns)
- Segmented copy number calls with confidence scores
- Primary output for CNV interpretation

### Visualization Plots
- Scatter plots showing copy number across chromosomes
- Segment overlays indicating called CNV regions
- PDF format for publication-quality figures

## Citation

If you use this module in your research, please cite:

**CNVkit:**
Talevich, E., Shain, A.H., Botton, T., & Bastian, B.C. (2016). CNVkit: Genome-Wide Copy Number Detection and Visualization from Targeted DNA Sequencing. *PLOS Computational Biology*, 12(4), e1004873.

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

- **[CNVkit Documentation](https://cnvkit.readthedocs.io/)**: Official CNVkit documentation
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
